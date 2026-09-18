"""Offline preflight and verified reuse of 19 unchanged cases for each model."""
from pathlib import Path
import shutil
import direct_raw_runner as r
from input_identity import sha256_file

RUNS={
    'luna':('direct_raw_41_v1','direct_raw_41_luna_idfixed_20260906'),
    'terra':('direct_raw_41_terra_v1','direct_raw_41_terra_idfixed_20260906'),
    'sol':('direct_raw_41_sol_v1','direct_raw_41_sol_idfixed_20260906'),
}
EXPECTED_CHANGED={'P8','P9','P10','P12','P13','P14','P15','P21','P22','P23','P24','P25','P28','P29','P30','P32','P33','P34','P35','P36','P37','P39'}

def main():
    _,specs=r.load_manifest(r.DEFAULT_MANIFEST)
    prompt_hash=r.sha256_text(r.DEFAULT_PROMPT.read_text(encoding='utf-8-sig').strip())
    audits={s.patient_id:r.build_patient_payload(s)[1] for s in specs}
    if any(a['remaining_exact_date_hits'] or not a['index_date_available'] for a in audits.values()): raise ValueError('Preflight failed')
    pending=[]
    for short,(old_name,new_name) in RUNS.items():
        old_dir=r.DEFAULT_OUTPUT_ROOT/old_name
        new_dir=r.DEFAULT_OUTPUT_ROOT/new_name
        if new_dir.exists(): raise ValueError(f'New output already exists: {new_dir}')
        reused=[]; changed=[]
        for spec in specs:
            source=old_dir/f'{spec.patient_id}_direct_raw.json'
            record=r.load_json(source)
            audit=audits[spec.patient_id]
            if record['prompt_sha256']!=prompt_hash or record['model_requested']!=f'gpt-5.6-{short}' or record['reasoning_effort']!='medium':
                raise ValueError('Original prompt/model/effort changed')
            if record['payload_sha256']==audit['payload_sha256']:
                r.validate_cached_output(record,spec,audit,prompt_hash,f'gpt-5.6-{short}','medium')
                reused.append({'patient_id':spec.patient_id,'source_file':str(source.resolve()),'source_sha256':sha256_file(source),'payload_sha256':audit['payload_sha256']})
            else: changed.append(spec.patient_id)
        if set(changed)!=EXPECTED_CHANGED or len(reused)!=19: raise ValueError(f'Unexpected change scope for {short}: {changed}')
        pending.append((short,new_dir,reused,changed))
    # No writes until all three preflights pass.
    for short,new_dir,reused,changed in pending:
        new_dir.mkdir()
        for row in reused:
            shutil.copyfile(row['source_file'],new_dir/f"{row['patient_id']}_direct_raw.json")
        r.write_json(new_dir/'reuse_provenance.json',{'status':'validated_reuse','reason':'Payload, prompt, model and reasoning unchanged; retain identical original output bytes.','reused':reused,'must_rerun':changed})
        snapshots=new_dir/'source_snapshots'
        snapshots.mkdir()
        for p in [r.DEFAULT_MANIFEST,r.DEFAULT_PROMPT,r.SCRIPT_DIR/'input_mapping_v2.json',r.SCRIPT_DIR/'direct_raw_runner.py',r.SCRIPT_DIR/'input_identity.py']:
            shutil.copyfile(p,snapshots/p.name)
        print(f'{short}: {len(reused)} reused, {len(changed)} to rerun -> {new_dir}',flush=True)
    r.write_json(r.SCRIPT_DIR/'audit/input_audit_idfixed_20260906.json',r.build_audit(r.DEFAULT_MANIFEST))

if __name__=='__main__': main()
