"""Synthetic upstream-merge -> OBER adapter -> R5 finalizer demonstration.

The article/clinical signals below are invented fixtures, not model outputs.
This validates the handoff contract, not an end-to-end clinical inference.
"""
from pathlib import Path
import argparse, hashlib, json, sys

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts/offline_guard'), str(ROOT/'RAG_re'), str(ROOT/'OBER_patient_delivery')]
import sitecustomize  # Explicitly activate the offline guard before research imports.
from rag_re.input_adapter import load_bundle, patient_id_from, select_candidate_pool, candidate_name, baseline_picked_names
from generate_r5_final_decisions import compose_decision

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'local_outputs/synthetic_handoff.json')
    args=parser.parse_args()
    source=ROOT/'examples/ober_handoff/merged.synthetic.json'
    root,bundle=load_bundle(source)
    pool,metadata=select_candidate_pool(root,bundle,{'merge_review_tiers':['review_high_priority','review_context_needed']})
    assert [candidate_name(c) for c in pool]==['Synthetic organism Alpha','Synthetic organism Beta']
    signal={
        'candidate_organism':'Synthetic organism Beta','canonical_organism':'Synthetic organism Beta',
        'review_tier':'review_high_priority','A1_strong':2,'A1_partial':0,'A1_match':2,
        'B_strict':False,'B_absolute_axis':True,'B_relative_axis':False,'B_site_aligned':True,
        'C_grade':'C0','D_state':'D_PASS','clinical_features_source':'invented synthetic fixture',
        'casefit_file':'invented fixture (no search, no LLM)','casefit_sha256':None,
        'latest_merge_file':'examples/ober_handoff/merged.synthetic.json',
        'latest_merge_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
    }
    result=compose_decision(cohort='synthetic',patient_id=patient_id_from(source,root,bundle),candidate_names=[candidate_name(c) for c in bundle['pathogen_candidates']],picked=[{'canonical_name':n,'rank':i+1} for i,n in enumerate(baseline_picked_names(bundle))],signals=[signal],source_snapshot={'synthetic':True,'merge_sha256':signal['latest_merge_sha256']})
    assert result['counts']['R5_rescued']==1
    assert [r['organism'] for r in result['selected']]==['Synthetic organism Alpha','Synthetic organism Beta']
    assert result['validation']['gold_answer_used_in_this_run'] is False
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8') as f:
        json.dump({'synthetic_only':True,'adapter_metadata':metadata,'decision':result},f,ensure_ascii=False,indent=2)
        f.write('\n')
    print('PASS: synthetic merge -> adapter -> R5, 1 preserved + 1 rescued; no network/API')
    return 0

if __name__=='__main__':raise SystemExit(main())
