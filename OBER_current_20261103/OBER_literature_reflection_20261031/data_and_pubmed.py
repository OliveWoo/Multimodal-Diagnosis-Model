"""Source-only preparation and organism-specific bidirectional PubMed retrieval."""
import argparse
import collections
import concurrent.futures
import json
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;RUN=HERE/'runs/reflection_v1'
sys.path.insert(0,str(ROOT/'OBER_clinical_framework_20261031'))
import experiment as clinical
import clinical_rules as names
from inspect_release import load,digest,require,write_new

SOURCE=ROOT/'OBER_clinical_framework_20261031/runs/full49_v1_verified_spans'
CUTOFF='2026/10/31'
SEARCH_NAMES={
    'Human alphaherpesvirus 1':['herpes simplex virus type 1','HSV-1'],
    'Human alphaherpesvirus 3':['varicella-zoster virus','VZV'],
    'Human betaherpesvirus 5':['cytomegalovirus','CMV'],
    'Human betaherpesvirus 6B':['human herpesvirus 6B','HHV-6B'],
    'Human betaherpesvirus 7':['human herpesvirus 7','HHV-7'],
    'Betapolyomavirus hominis':['BK polyomavirus'],
    'Severe acute respiratory syndrome-related coronavirus':['SARS-CoV-2'],
    'Nakaseomyces glabratus':['Candida glabrata'],
}


def save(path,data):write_new(path,clinical.sources.old.encoded(data))


def prepare():
    clinical.RUN=SOURCE;source_plan=clinical.verify_outputs()
    raw_rows=clinical.sources.v4.all_rows()
    components={(r['case_id'],r['organism']):r for r in load(SOURCE/'candidates_not_model_input.json')}
    patients={p['case_id']:p for p in source_plan['patients']};candidates=[]
    for row in raw_rows:
        rec,c=row['record'],row['candidate']
        if c['reporting']['tier']!='Context':continue
        cid=rec['case_id'];p=patients[cid];common=load(SOURCE/p['folder']/'input.json');ledger=load(SOURCE/p['folder']/'ledger.json')
        by_id={f['id']:f for f in common['facts']}
        relevant=names.matched_observations(c['organism_name'],ledger)
        # Every recorded positive alternative remains visible; no ranking against
        # upstream picks or reference. Retain candidate negatives/equivocal too.
        obs_ids={o['fact_id'] for o in relevant}|{o['fact_id'] for o in ledger['observations'] if o['detection']=='positive'}
        facts=[f for f in common['facts'] if f['kind'] in ['index','background','imaging'] or f['id'] in obs_ids]
        token=digest((cid+'|'+c['organism_key']).encode())[:14]
        component=components[(cid,c['organism_name'])]
        typed={b['seq_id']:b.get('nucleic_type') for b in c['analytical_evidence'].get('per_test_signals') or []}
        for i,b in enumerate(component['branches'],1):
            value={**b,'nucleic_type':typed[b['seq_id']],
                'interpretation':'Analytical measurements, not calibrated causality thresholds; technical branches are not independent confirmation.'}
            facts.append({'id':f'M_{token}_{i:02}','kind':'mngs','temporal_role':b['temporal_role'],
                'text':clinical.flatten(value)})
        premises=[]
        for field,seq in [('reporting/blockers',c['reporting'].get('blockers') or []),
                          ('reporting/cautions',c['reporting'].get('cautions') or []),
                          ('reporting_packet/reasons',row['detail'].get('reasons') or [])]:
            for i,s in enumerate(seq):
                if s and s not in [v['text'] for v in premises]:
                    premises.append({'id':f'U{len(premises)+1:02}','text':s,'source_field':field+'/'+str(i)})
        if not premises:premises=[{'id':'U00','text':'Retained as Context; no explicit reason supplied. Do not invent one.','source_field':'reporting/tier'}]
        folder=RUN/'cases'/cid/token
        packet={'candidate_token':token,'organism':c['organism_name'],
            'upstream_premises':premises,'facts':facts,
            'source_limits':{'available_source_facts':len(common['facts']),'provided_source_facts':len(facts)-len(component['branches']),
                'withheld_microbiology':'Noncandidate negatives and nonresult records; not proof of absence.',
                'time':'index +/-2 calendar dates; unknown-date context is not contemporaneous confirmation.',
                'sampling':'No validated reads/RPM threshold; no complete EMR; original positive decisions not supplied.'}}
        require(not any(clinical.sources.BLOCK.search(f['text']) for f in facts),'Unsafe clinical source')
        save(folder/'packet.json',packet)
        save(folder/'source_link_not_model_input.json',{'case_id':cid,'original_case_folder':str((SOURCE/p['folder']).relative_to(ROOT)),
            'source_provenance_sha256':digest((SOURCE/p['folder']/'source_provenance.json').read_bytes()),
            'upstream_result':rec['patient_result_file'],'upstream_packet':rec['reporting_packet_file'],
            'upstream_result_sha256':rec['patient_result_sha256'],'upstream_packet_sha256':rec['reporting_packet_sha256'],
            'new_mngs_component_source':str((SOURCE/'candidates_not_model_input.json').relative_to(ROOT)),
            'organism_key':c['organism_key']})
        candidates.append({'case_id':cid,'group':rec['analysis_group'],'organism':c['organism_name'],
            'organism_key':c['organism_key'],'token':token,'folder':folder.relative_to(RUN).as_posix()})
    require(len(candidates)==319 and len({r['case_id'] for r in candidates})==49,'Cohort mismatch')
    organisms=sorted({r['organism'] for r in candidates})
    save(RUN/'source_plan.json',{'cases':49,'candidates':candidates,'organisms':organisms,
        'publication_cutoff':CUTOFF,'baseline_positives_preserved':97,'reference_read':False,
        'source_output_freeze':digest((SOURCE/'output_freeze.json').read_bytes()),
        'data_adapter_sha256':digest(Path(__file__).read_bytes())})
    print(json.dumps({'prepared_cases':49,'candidates':len(candidates),'unique_organisms':len(organisms)},ensure_ascii=False))


_lock=threading.Lock();_last=0.0
def get(url):
    global _last
    for attempt in range(3):
        with _lock:
            delay=max(0,.38-(time.monotonic()-_last))
            if delay:time.sleep(delay)
            _last=time.monotonic()
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'OBER-literature-reflection-research/1.0'})
            with urllib.request.urlopen(req,timeout=50) as response:return response.read()
        except (urllib.error.URLError,TimeoutError):
            if attempt==2:raise
            time.sleep(1+attempt)


def queries(name):
    terms=[name]+SEARCH_NAMES.get(name,[])
    organism='('+' OR '.join('"'+s.replace('"',' ')+'"[Title/Abstract]' for s in terms)+')'
    base=organism+' AND (pneumonia[Title/Abstract] OR pulmonary[Title/Abstract] OR respiratory[Title/Abstract])'
    date=f' AND (1900/01/01:{CUTOFF}[dp])'
    return {'general':base+date,'counter':base+' AND (colonization[Title/Abstract] OR colonisation[Title/Abstract] OR contamination[Title/Abstract] OR reactivation[Title/Abstract] OR shedding[Title/Abstract] OR asymptomatic[Title/Abstract])'+date}


def retrieve_one(name):
    folder=RUN/'literature'/digest(name.encode())[:16]
    if (folder/'pubmed.json').exists():return load(folder/'pubmed.json')
    searches={};base='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    for lane,query in queries(name).items():
        path=folder/(lane+'_search.json')
        recovery=folder/(lane+'_search_recovered.json')
        result=load(recovery) if recovery.exists() else load(path) if path.exists() else None
        if result is None or 'count' not in result.get('result',{}):
            for attempt in range(3):
                data=get(base+'esearch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','term':query,
                    'retmode':'json','retmax':10,'sort':'relevance','tool':'ober_literature_reflection'}))
                current={'query':query,'result':json.loads(data)['esearchresult']}
                if 'count' in current['result']:
                    result=current;save(recovery if path.exists() else path,result);break
                save(folder/(lane+'_backend_error_'+digest(data)[:12]+'.json'),current)
                time.sleep(1+attempt)
            require(result is not None and 'count' in result.get('result',{}),'PubMed search backend failed')
        require(result['query']==query,'Changed retrieval cache');searches[lane]=result
    general=searches['general']['result'].get('idlist') or [];counter=searches['counter']['result'].get('idlist') or []
    ids=list(general[:8]);counter_added=[]
    for pmid in counter:
        if pmid not in ids and len(counter_added)<2:ids.append(pmid);counter_added.append(pmid)
    for pmid in general+counter:
        if len(ids)<10 and pmid not in ids:ids.append(pmid)
    articles=[]
    if ids:
        xml=folder/'pubmed_raw.xml'
        if xml.exists():raw=xml.read_bytes()
        else:
            raw=get(base+'efetch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','id':','.join(ids),'retmode':'xml','tool':'ober_literature_reflection'}));write_new(xml,raw)
        def txt(x):return ' '.join(''.join(x.itertext()).split()) if x is not None else ''
        def parse(data):
            by_id={}
            root=ET.fromstring(data)
            for a in list(root):
                if a.tag not in ['PubmedArticle','PubmedBookArticle']:continue
                pmid=txt(a.find('.//MedlineCitation/PMID')) or txt(a.find('.//BookDocument/PMID'))
                by_id[pmid]={'pmid':pmid,'title':txt(a.find('.//ArticleTitle')) or txt(a.find('.//BookTitle')),
                    'abstract':' '.join(txt(x) for x in a.findall('.//Abstract/AbstractText')),
                    'publication_types':[txt(x) for x in a.findall('.//PublicationTypeList/PublicationType')],
                    'record_kind':a.tag,'url':'https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/',
                    'counter_lane_selected':pmid in counter_added}
            return by_id
        by_id=parse(raw)
        if not set(ids)<=set(by_id):
            # A transient failed search may have changed the selected ID list.
            # Preserve the first XML; fetch the explicitly selected replacement.
            fresh=folder/('pubmed_raw_recovered_'+digest(','.join(ids).encode())[:12]+'.xml')
            if fresh.exists():raw=fresh.read_bytes()
            else:
                raw=get(base+'efetch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','id':','.join(ids),'retmode':'xml','tool':'ober_literature_reflection'}));write_new(fresh,raw)
            by_id=parse(raw)
        require(set(ids)<=set(by_id),'PubMed XML missing selected article')
        articles=clinical.sources.old.with_sentences([by_id[i] for i in ids])
    out={'organism':name,'queries':queries(name),'articles':articles,'selected_pmids':ids,
        'general_hits':int(searches['general']['result']['count']),
        'counter_hits':int(searches['counter']['result']['count']),'technical_success':True,
        'publication_cutoff':CUTOFF,'max_articles':10,'abstract_count':sum(bool(a['abstract']) for a in articles)}
    save(folder/'pubmed.json',out)
    print(json.dumps({'retrieved':name,'articles':len(articles),'abstracts':out['abstract_count']},ensure_ascii=True),flush=True)
    return out


def retrieve():
    plan=load(RUN/'source_plan.json');errors=[];results={}
    def work(name):
        try:results[name]=retrieve_one(name)
        except Exception as exc:
            errors.append({'organism':name,'type':type(exc).__name__,'error':str(exc)[:180]});print(json.dumps({'retrieval_failed':name,'type':type(exc).__name__}),flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(work,plan['organisms']))
    require(not errors,'Retrieval technical failures, not zero articles: '+json.dumps(errors))
    save(RUN/'retrieval_summary.json',{'organisms':len(results),'article_pairs':sum(len(r['articles']) for r in results.values()),
        'with_abstract':sum(r['abstract_count']>0 for r in results.values()),'technical_failures':0,
        'per_organism':{n:{k:r[k] for k in ['selected_pmids','general_hits','counter_hits','abstract_count']} for n,r in results.items()}})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','retrieve']);{'prepare':prepare,'retrieve':retrieve}[p.parse_args().command]()
