"""Three broad, answer-independent checks on frozen literature reflections."""
import itertools
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'OBER_reason_reflection_20261031'))
import reason_common as old
import reason_judge as judge
import native_signal as native

RUN=HERE/'runs/modules_v1'
PARENT=ROOT/'OBER_reason_reflection_20261031/runs/reason_v1'
load,save,require,digest,write_new=old.load,old.save,old.require,old.digest,old.write_new
PRIMARY='M_C1_Q0_S0_BOTH'


def freeze(stage,paths):
    save(RUN/(stage+'_freeze.json'),{'files':{str(p.resolve()):digest(p.read_bytes()) for p in paths},
        'posthoc_development_not_independent_test':True})


def verify(stage):
    for name,h in load(RUN/(stage+'_freeze.json'))['files'].items():
        require(digest(Path(name).read_bytes())==h,'Changed frozen module file '+name)


def parent_verify(stage):
    for name,h in load(PARENT/(stage+'_freeze.json'))['files'].items():
        require(digest(Path(name).read_bytes())==h,'Changed frozen parent file '+name)


def settings():
    return [{'id':f'M_C{cc}_Q{qq}_S{ss}_{signal.upper()}',
        'C':bool(cc),'Q':bool(qq),'S':bool(ss),'signal':signal}
        for cc,qq,ss in itertools.product([0,1],repeat=3) for signal in ['none','either','both']]


def feature(q,p,gq,gp):
    base=judge.selection(q,'broad') if q else {'add':False,'qualifying_pmids':[]}
    if not q:return {'T_broad':False,'G_broad':bool(gq and judge.selection(gq,'broad')['add']),
        'C':False,'Q':False,'S':False,'unavailable_T':True,'closure_pmids':[],
        'shared_clinical_fact_ids':[],'counter_severity':None,'papers_same':p['papers']==gp['papers']}
    def cited_abstract(x):
        papers={a['pmid']:a for a in p['papers']}
        return {r['pmid'] for r in x['paper_refs'] if r['pmid'] in papers
            and judge.has_abstract_citation(papers[r['pmid']],r['sentence_ids'])}
    good=set(base['qualifying_pmids'])
    closure=good & cited_abstract(q['case_link']) & cited_abstract(q['revision'])
    facts={f['id']:f for f in p['patient_facts']}
    common=set(q['case_link']['patient_fact_ids']) & set(q['revision']['patient_fact_ids'])
    clinical=sorted(i for i in common if facts[i]['kind'] in ['background','microbiology','imaging'])
    gadd=bool(gq and judge.selection(gq,'broad')['add'])
    return {'T_broad':base['add'],'G_broad':gadd,
        'C':q['counterevidence']['severity'] in ['none','limited_uncertainty'],
        'Q':bool(closure and clinical),'S':gadd,'unavailable_T':False,
        'closure_pmids':sorted(closure),'shared_clinical_fact_ids':clinical,
        'counter_severity':q['counterevidence']['severity'],'papers_same':p['papers']==gp['papers']}


def select(f,flags,spec):
    signal=spec['signal']=='none' or native.passes(flags,spec['signal'])
    return bool(f['T_broad'] and signal and all(not spec[k] or f[k] for k in ['C','Q','S']))


def prepare():
    parent_verify('judge_inputs');parent_verify('judge_outputs')
    require(not (RUN/'plan.json').exists(),'Prepared run already exists')
    for name,h in load(native.NRUN/'prediction_freeze.json')['files'].items():
        require(digest((native.NRUN/name).read_bytes())==h,'Native source changed')
    source=load(PARENT/'source_plan.json')['candidates']
    features={x['token']:x for x in load(native.NRUN/'native_features_not_model_input.json')}
    rows=[]
    for row in source:
        inputs={a:load(PARENT/'judgements'/row['case_id']/row['token']/a/'input.json') for a in ['G','T']}
        outputs={}
        for arm,p in inputs.items():
            folder=PARENT/'judgements'/row['case_id']/row['token']/arm
            outputs[arm]=old.validate_saved(folder,judge.validate) if any(a['abstract'] for a in p['papers']) else None
        f=feature(outputs['T'],inputs['T'],outputs['G'],inputs['G'])
        flags=features[row['token']]
        rows.append({k:row[k] for k in ['case_id','group','token','organism']}|{'features':f,
            'native_flags':{k:flags[k] for k in ['absolute','relative','native_rpm_threshold','native_rank_band']},
            'decisions':{s['id']:select(f,flags,s) for s in settings()}})
    save(RUN/'predictions_without_answers.json',rows)
    save(RUN/'plan.json',{'primary':PRIMARY,'settings':settings(),'candidate_count':len(rows),
        'posthoc_development_not_independent_test':True,'no_new_API_calls':True,'no_new_patient_information':True,
        'no_gold_loaded_during_prepare':True,'historical_labels_previously_seen':True,
        'parent_judge_input_freeze_sha256':digest((PARENT/'judge_inputs_freeze.json').read_bytes()),
        'parent_judge_output_freeze_sha256':digest((PARENT/'judge_outputs_freeze.json').read_bytes()),
        'native_prediction_freeze_sha256':digest((native.NRUN/'prediction_freeze.json').read_bytes())})
    freeze('prediction',[RUN/'predictions_without_answers.json',RUN/'plan.json',Path(__file__),HERE/'PROTOCOL_zh.md',
        ROOT/'OBER_reason_reflection_20261031/reason_judge.py',ROOT/'OBER_reason_reflection_20261031/native_signal.py'])
    print({'frozen_candidates':len(rows),'fixed_settings':len(settings()),'primary':PRIMARY,'new_API_calls':0},flush=True)


if __name__=='__main__':prepare()
