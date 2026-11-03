"""PubMed question-axis retrieval; only controlled public concepts leave locally."""
import argparse
import collections
import concurrent.futures
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
import reason_common as c

BASE='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
AXIS_TERMS={
    'causality':'(colonization[tiab] OR colonisation[tiab] OR pathogenicity[tiab] OR causative[tiab] OR "clinical significance"[tiab] OR "diagnostic criteria"[tiab])',
    'specimen_localization':'(blood[tiab] OR plasma[tiab] OR bronchoalveolar[tiab] OR sputum[tiab] OR "cell-free DNA"[tiab] OR dissemination[tiab])',
    'host_syndrome':'("risk factors"[tiab] OR host[tiab] OR immunocompromised[tiab] OR "clinical presentation"[tiab] OR diagnosis[tiab])',
    'test_interpretation':'(metagenomic[tiab] OR mNGS[tiab] OR "diagnostic accuracy"[tiab] OR "diagnostic criteria"[tiab] OR "molecular diagnosis"[tiab] OR PCR[tiab])'}
HOST_TERMS={
    'immunocompromised':'(immunocompromised[tiab] OR immunosuppressed[tiab] OR corticosteroid*[tiab] OR transplant*[tiab] OR hematologic*[tiab])',
    'structural_airway':'(bronchiectasis[tiab] OR "cystic fibrosis"[tiab] OR COPD[tiab] OR "chronic lung"[tiab])',
    'critical_illness_only':'("critical illness"[tiab] OR ICU[tiab] OR ventilat*[tiab] OR "intensive care"[tiab])',
    'other_or_general':'','unknown':''}
COUNTER='(colonization[tiab] OR colonisation[tiab] OR contamination[tiab] OR reactivation[tiab] OR asymptomatic[tiab] OR "false positive"[tiab] OR bystander[tiab])'
_request_lock=threading.Lock()
_last_request=0.0


def pubmed_get(url):
    global _last_request
    for attempt in range(5):
        with _request_lock:
            delay=max(0,1.1-(time.monotonic()-_last_request))
            if delay:time.sleep(delay)
            _last_request=time.monotonic()
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'OBER-reason-reflection-research/1.0'})
            with urllib.request.urlopen(req,timeout=50) as response:return response.read()
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as exc:
            status=getattr(exc,'code',None)
            if attempt==4 or status is not None and status not in [408,429,500,502,503,504,520,522,524]:raise
            c.save(c.RUN/'pubmed_transport'/f'error_{time.time_ns()}.json',{
                'error_type':type(exc).__name__,'status':status,'attempt':attempt+1,
                'url_sha256':c.digest(url.encode()),'not_no_evidence':True})
            time.sleep(min(40,5*2**attempt))


def clean_article(a):
    result={k:a[k] for k in ['pmid','title','abstract','sentences']}
    return result|{'publication_types':a.get('publication_types') or [],
        'record_kind':a.get('record_kind','not_recorded_in_original_cache'),
        'url':a.get('url') or 'https://pubmed.ncbi.nlm.nih.gov/'+a['pmid']+'/'}


def queries(organism,agenda):
    # Nomenclature aliases only: never import retrospective scoring name mappings.
    names=[organism]+c.old.d.SEARCH_NAMES.get(organism,[])
    term='('+' OR '.join('"'+n.replace('"',' ')+'"[tiab]' for n in names)+')'
    base=term+' AND (pneumonia[tiab] OR pulmonary[tiab] OR respiratory[tiab])'
    date=' AND (1900/01/01:2026/10/31[dp]) AND hasabstract[text]'
    result={}
    for qid,q in agenda['questions'].items():
        core=base+' AND '+AXIS_TERMS[q['axis']]
        # The actual case population is a retrieval context, not a selection rule.
        # Keep an axis-only fallback when this controlled context returns no hits.
        host=HOST_TERMS[agenda['host_context']['category']]
        result[qid]={'primary':core+(' AND '+host if host else '')+date,
            'fallback':core+date,'axis':q['axis']}
    result['counter']={'primary':base+' AND '+COUNTER+date,
        'fallback':base+' AND '+COUNTER+date,'axis':'counter'}
    return result


def search(query):
    folder=c.RUN/'searches'/c.digest(query.encode())[:20];path=folder/'search.json'
    if path.exists():
        result=c.load(path);c.require(result['query']==query,'Query cache collision');return result
    for attempt in range(3):
        data=pubmed_get(BASE+'esearch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','term':query,
            'retmode':'json','retmax':8,'sort':'relevance','tool':'ober_reason_reflection'}))
        import json
        body=json.loads(data).get('esearchresult',{})
        if str(body.get('count','')).isdigit():
            result={'query':query,'result':body,'sent_patient_data':False};c.save(path,result);return result
        c.save(folder/f'backend_error_{time.time_ns()}.json',body)
        time.sleep(attempt+1)
    raise ValueError('PubMed backend failure, not zero articles')


def select_ids(general,searches):
    ids=[];lanes=collections.defaultdict(list)
    def add(seq,lane,limit):
        added=0
        for pmid in seq:
            if lane not in lanes[pmid]:lanes[pmid].append(lane)
            if pmid not in ids and added<limit and len(ids)<10:ids.append(pmid);added+=1
    add(general,'general_anchor',2)
    add(searches['Q1'],'Q1',3);add(searches['Q2'],'Q2',3);add(searches['counter'],'counter',2)
    for lane in ['Q1','Q2','counter']:add(searches[lane],lane,10)
    add(general,'general_fill',10)
    return ids,{p:lanes[p] for p in ids}


def parse_articles(raw):
    def txt(x):return ' '.join(''.join(x.itertext()).split()) if x is not None else ''
    articles=[]
    for a in ET.fromstring(raw):
        if a.tag not in ['PubmedArticle','PubmedBookArticle']:continue
        pmid=txt(a.find('.//MedlineCitation/PMID')) or txt(a.find('.//BookDocument/PMID'))
        articles.append({'pmid':pmid,'title':txt(a.find('.//ArticleTitle')) or txt(a.find('.//BookTitle')),
            'abstract':' '.join(txt(x) for x in a.findall('.//Abstract/AbstractText')),
            'publication_types':[txt(x) for x in a.findall('.//PublicationTypeList/PublicationType')],
            'record_kind':a.tag,'url':'https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/'})
    return c.old.d.clinical.sources.old.with_sentences(articles)


def retrieve():
    c.verify('planning_inputs');c.verify('planning_outputs')
    if (c.RUN/'retrieval_outputs_freeze.json').exists():c.verify('retrieval_outputs');print('Retrieval complete; no requests.');return
    source=c.load(c.RUN/'source_plan.json');lookup={};jobs={}
    for row in source['candidates']:
        agenda=c.load(c.RUN/row['folder']/'result.json')['output'];qs=queries(row['organism'],agenda)
        lookup[row['token']]=qs
        for q in qs.values():jobs[q['primary']]=q
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        result=list(pool.map(search,jobs))
    indexed={r['query']:r for r in result};fallbacks={}
    for qs in lookup.values():
        for q in qs.values():
            if int(indexed[q['primary']]['result']['count'])==0 and q['fallback']!=q['primary']:
                fallbacks[q['fallback']]=True
    if fallbacks:
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            for r in pool.map(search,fallbacks):indexed[r['query']]=r
    all_ids=set();by_token={};old_articles={}
    for row in source['candidates']:
        lib=c.load(c.SOURCE/'literature'/c.digest(row['organism'].encode())[:16]/'pubmed.json')
        general=lib['articles'];found={};history={}
        for a in general:old_articles.setdefault(a['pmid'],clean_article(a))
        for lane,q in lookup[row['token']].items():
            chosen=q['primary']
            if int(indexed[chosen]['result']['count'])==0 and q['fallback']!=chosen:chosen=q['fallback']
            found[lane]=indexed[chosen]['result'].get('idlist') or []
            history[lane]={'axis':q['axis'],'used_query':chosen,'hit_count':int(indexed[chosen]['result']['count'])}
        selected,lanes=select_ids([a['pmid'] for a in general],found)
        by_token[row['token']]={'G':general,'selected':selected,'query_history':history,'lanes':lanes}
        all_ids.update(selected)
    cache=c.RUN/'article_cache'
    for pmid,a in old_articles.items():c.save(cache/(pmid+'.json'),a)
    missing=sorted(p for p in all_ids if not (cache/(p+'.json')).exists())
    chunks=[missing[start:start+40] for start in range(0,len(missing),40)]
    def fetch(ids):
        key=c.digest(','.join(ids).encode())[:20];parsed={}
        for attempt in range(3):
            suffix='' if attempt==0 else '_recovery'+str(attempt)
            path=c.RUN/'efetch'/(key+suffix+'.xml')
            if path.exists():raw=path.read_bytes()
            else:
                raw=pubmed_get(BASE+'efetch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','id':','.join(ids),
                    'retmode':'xml','tool':'ober_reason_reflection'}));c.write_new(path,raw)
            try:parsed={a['pmid']:a for a in parse_articles(raw)}
            except ET.ParseError:parsed={}
            if set(ids)<=set(parsed):break
        c.require(set(ids)<=set(parsed),'EFetch missing selected records, not zero evidence')
        for pmid in ids:c.save(cache/(pmid+'.json'),clean_article(parsed[pmid]))
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(fetch,chunks))
    files=[];summary=collections.Counter()
    for row in source['candidates']:
        lib=by_token[row['token']];targeted=[c.load(cache/(pmid+'.json')) for pmid in lib['selected']]
        path=c.RUN/'retrieval'/row['case_id']/row['token']/'papers.json'
        c.save(path,{'G':lib['G'],'T':targeted,'query_history':lib['query_history'],
            'selected_lanes_not_stance':lib['lanes'],'max_per_candidate_per_arm':10,
            'pubmed_received_only_controlled_terms':True})
        files.append(path)
        for arm,papers in [('G',lib['G']),('T',targeted)]:
            summary[arm+'_article_pairs']+=len(papers)
            summary[arm+'_with_abstract']+=any(a['abstract'] for a in papers)
        summary['both_with_abstract']+=any(a['abstract'] for a in lib['G']) and any(a['abstract'] for a in targeted)
    files+=list((c.RUN/'searches').rglob('search.json'))+list((c.RUN/'efetch').glob('*.xml'))+list(cache.glob('*.json'))
    c.save(c.RUN/'retrieval_summary.json',dict(summary)|{'unique_searches':len(indexed),
        'unique_targeted_pmids':len(all_ids),'newly_fetched_pmids':len(missing),'technical_failures':0})
    c.freeze('retrieval_outputs',files+[c.RUN/'retrieval_summary.json']+c.deps('reason_retrieval.py'))
    print(dict(summary)|{'unique_searches':len(indexed),'newly_fetched_pmids':len(missing),'no_labels_read':True})


if __name__=='__main__':retrieve()
