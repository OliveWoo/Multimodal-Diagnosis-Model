"""Two fresh QC consensus blocks, with explicitly shared fixed clinical controls."""
import argparse
import copy
import json
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'OBER_rule_QC_20261103'))
import qc_common as s
import qc_runner as runner
import qc_review as review
import qc_rules as rules
import qc_eval

PARENT=ROOT/'OBER_rule_QC_20261103/runs/qc_v1'
BUNDLE=HERE/'runs/repeat_v1'
ROUNDS={x:BUNDLE/('round_'+x) for x in ['A','B']}
COMMON=BUNDLE/'common_controls'


def verify_file_manifest(p):
    for name,h in s.load(p)['files'].items():s.require(s.digest(Path(name).read_bytes())==h,'Changed parent/control file '+name)


def freeze_bundle(stage,paths):
    s.save(BUNDLE/(stage+'_freeze.json'),{'files':{str(p.resolve()):s.digest(p.read_bytes()) for p in sorted(set(paths))},
        'two_fresh_QC_blocks_not_full_OBER_reruns':True,'independent_test':False})


def bind(label):
    folder=COMMON if label=='common' else ROUNDS[label]
    s.RUN=folder;s.c.RUN=folder;s.HERE=folder
    return folder


def prepare():
    global_path=BUNDLE/'setup_freeze.json'
    if global_path.exists():verify_file_manifest(global_path);print('Repeat setup already frozen.');return
    for stage in ['inputs','QC_outputs','article_decisions','A1_inputs','A1_outputs','predictions']:
        verify_file_manifest(PARENT/(stage+'_freeze.json'))
    parent=s.load(PARENT/'plan.json');paths=[]
    for label in ['A','B']:
        f=bind(label)
        s.require(not list(f.glob('articles/*/r*/result.json')) and not list(f.glob('articles/*/r*/*_raw_response.json')),'Fresh test directory contains prior QC results')
        plan=copy.deepcopy(parent);plan['fresh_repeat_block']=label;plan['G_T_and_A1_fixed_controls']=True
        for a in plan['articles']:
            for name in ['source.json','routing.json']:
                dest=f/a['folder']/name;s.write_new(dest,(PARENT/a['folder']/name).read_bytes());paths.append(dest)
            packet=s.load(f/a['folder']/'source.json')
            s.require(rules.route(packet)==s.load(f/a['folder']/'routing.json'),'Routing changed')
        for j in plan['QC_jobs']:
            for name in ['input.json','schema.json']:
                dest=f/j['folder']/name;s.write_new(dest,(PARENT/j['folder']/name).read_bytes());paths.append(dest)
        for row in plan['potential']:
            for name in ['a1_patient_view.json','papers.json']:
                dest=f/row['folder']/name;s.write_new(dest,(PARENT/row['folder']/name).read_bytes());paths.append(dest)
        s.save(f/'plan.json',plan);paths.append(f/'plan.json')
        deps=list(HERE.glob('*.py'))+[HERE/'PROTOCOL_zh.md',HERE/'README.md']
        deps += [Path(m.__file__) for m in [s,runner,review,rules,qc_eval,s.c]]
        s.freeze('inputs',[p for p in paths if f in p.parents]+deps)
        paths.append(f/'inputs_freeze.json')
    s.require(len(parent['QC_jobs'])==114,'Unexpected repeat batch size')
    freeze_bundle('setup',paths+list(HERE.glob('*.py'))+[HERE/'PROTOCOL_zh.md',HERE/'README.md',PARENT/'plan.json',PARENT/'QC_outputs_freeze.json'])
    print({'new_QC_calls_planned':228,'fresh_rounds':['A','B'],'old_QC_outputs_copied':0,'answers_not_read':True},flush=True)


def prepare_common_A1():
    if (COMMON/'A1_inputs_freeze.json').exists():verify_file_manifest(COMMON/'A1_inputs_freeze.json');return
    parent=s.load(PARENT/'plan.json');jobs=copy.deepcopy(s.load(PARENT/'A1_jobs.json'));existing={j['token'] for j in jobs};needed=[]
    for row in parent['potential']:
        required=row['old_S']
        for label in ['A','B']:
            bind(label)
            required=required or any(all(g['pass'] for g in runner.article_gate(row,strategy).values()) for strategy in runner.STRATEGIES)
        if required and row['token'] not in existing:needed.append(row)
    files=[];bind('common')
    for j in jobs:
        j['reused']=True;f=Path(j['folder']);files += [f/'input.json',f/'schema.json']
        files += s.response_files(f) if j['has_abstract'] else [f/'unavailable.json']
    for row in needed:
        view=s.load(PARENT/row['folder']/'a1_patient_view.json')
        for paper in s.load(PARENT/row['folder']/'papers.json'):
            for rep in range(2):
                f=COMMON/'alignments'/row['case_id']/row['token']/paper['pmid']/f'r{rep}'
                p={'organism':row['organism'],'case_card':view,'article':paper}
                aliases=s.a.legacy.CLINICAL_NAMES.get(row['organism'],[])
                if aliases:p['documented_clinical_names_not_scoring_aliases']=aliases
                s.save(f/'input.json',p);s.save(f/'schema.json',s.olda1.schema(p))
                jobs.append({k:row[k] for k in ['token','case_id','organism']}|{'pmid':paper['pmid'],'replicate':rep,
                    'folder':str(f.resolve()),'reused':False,'has_abstract':bool(paper['abstract'])})
                files += [f/'input.json',f/'schema.json']
    s.save(COMMON/'A1_jobs.json',jobs);s.freeze('A1_inputs',files+[COMMON/'A1_jobs.json'])
    print({'common_new_A1_calls_if_needed':sum(not j['reused'] and j['has_abstract'] for j in jobs),
        'same_fixed_A1_data_for_both_rounds':True},flush=True)


def attach_common_A1(label):
    f=bind(label)
    if (f/'A1_outputs_freeze.json').exists():s.verify('A1_inputs');s.verify('A1_outputs');return
    shared=s.load(COMMON/'A1_jobs.json');jobs=[j|{'reused':True} for j in shared]
    s.save(f/'A1_jobs.json',jobs)
    inputs=[f/'A1_jobs.json',COMMON/'A1_inputs_freeze.json',COMMON/'A1_outputs_freeze.json'];outputs=[]
    for j in jobs:
        d=Path(j['folder']);inputs += [d/'input.json',d/'schema.json']
        outputs += s.response_files(d) if j['has_abstract'] else [d/'unavailable.json']
    s.freeze('A1_inputs',inputs);s.freeze('A1_outputs',outputs)


def pipeline():
    prepare();verify_file_manifest(BUNDLE/'setup_freeze.json')
    for label in ['A','B']:
        bind(label);runner.run_QC();runner.derive_articles()
        print({'fresh_block_completed':label},flush=True)
    prepare_common_A1();bind('common');runner.run_A1()
    for label in ['A','B']:
        attach_common_A1(label);runner.predict()
    freeze_bundle('both_predictions',[ROUNDS[x]/'predictions_freeze.json' for x in ['A','B']]+[COMMON/'A1_outputs_freeze.json'])
    print('Both fresh QC blocks and all predictions frozen; no reference labels read.',flush=True)


def round_usage():
    import collections
    totals=collections.Counter();models=collections.Counter()
    for p in s.RUN.glob('articles/*/r*/*_raw_response.json'):
        q=s.load(p)['response'];totals['saved_responses']+=1;models[q.get('model','unknown')]+=1
        for k in ['input_tokens','output_tokens','total_tokens']:totals[k]+=(q.get('usage') or {}).get(k,0)
    out={'new_QC_response_usage':dict(totals),'actual_models':dict(models),'fixed_A1_controls_shared_between_rounds':True,
        'common_A1_control_manifest':str(COMMON/'A1_outputs_freeze.json'),'not_invoice':True}
    s.save(s.RUN/'api_usage.json',out);return out


def evaluate():
    verify_file_manifest(BUNDLE/'setup_freeze.json');verify_file_manifest(BUNDLE/'both_predictions_freeze.json')
    original=qc_eval.usage
    try:
        qc_eval.usage=round_usage
        for label in ['A','B']:bind(label);qc_eval.main()
    finally:qc_eval.usage=original
    from compare_runs import compare
    compare()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','pipeline','evaluate'])
    {'prepare':prepare,'pipeline':pipeline,'evaluate':evaluate}[p.parse_args().stage]()
