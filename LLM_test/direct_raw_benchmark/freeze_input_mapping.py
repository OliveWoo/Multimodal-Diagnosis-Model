"""Freeze the reviewed identity joins using raw sources only; no Gold or OBER reads."""
import hashlib
import json
from pathlib import Path

BASE=Path(__file__).resolve().parent
SECTIONS=('admission_diagnosis','underlying','CBC','other_lab','culture','filmarray','gm_test','molecular_microbiology','image')
# Explicit audited joins, NOT a numeric +/-1 heuristic. P16/P20 identify specific
# specimens in one-to-many case codes; retain those already-matching inputs.
KH_GLOBAL={3:3,4:4,6:6,8:7,9:8,10:9,12:11,13:12,14:13,15:14,16:16,17:17,18:18,19:19,20:20,21:22,22:23,23:24,24:25,25:26,28:29,29:30,30:31,32:33,33:34,34:35,35:36,36:37,37:38,39:40}
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def build_mapping():
    manifest=read(BASE/'cohort_manifest.json')
    index_file=BASE/'../data/patient_info_standardized/legacy_patient_index.json'
    index=read(index_file)
    by_id={int(r['legacy_patient_id'].rsplit('_',1)[1]):r for r in index}
    raw_root=(BASE/manifest['raw_mngs_source']).resolve()
    rows=[]
    for cohort in manifest['cohorts']:
        clinical_root=(BASE/cohort['clinical_root']).resolve()
        for patient in cohort['patient_ids']:
            pid=int(patient[1:])
            gid=KH_GLOBAL[pid] if cohort['name'].startswith('kmuh') else pid
            entry=by_id[gid]
            expected=f'K{pid:03}' if cohort['name'].startswith('kmuh') else entry['case_code']
            if entry['case_code']!=expected: raise ValueError(f'Unexpected case code for {patient}')
            raw=raw_root/f'NGS_patient_{gid}_all_RK_NTC_microbes.json'
            raw_data=read(raw)
            if len(raw_data)!=1 or raw_data[0]['specimen_code']!=entry['specimen_id']:
                raise ValueError(f'Expected one explicit indexed specimen for {patient}')
            clinical=clinical_root/f'NGS_patient_{pid}_json'
            hashes={s:sha(clinical/f'NGS_patient_{pid}_{s}.json') if (clinical/f'NGS_patient_{pid}_{s}.json').exists() else None for s in SECTIONS}
            notes=[]
            if pid in (16,20) and cohort['name'].startswith('kmuh'): notes.append('Case code has multiple index rows; this specimen is explicitly selected and agrees with the existing benchmark input.')
            if pid==32 and cohort['name'].startswith('kmuh'): notes.append('KH age/sex unavailable; no demographic data imputed.')
            if pid==33 and cohort['name'].startswith('kmuh'): notes.append('KH vs standardized age discrepancy remains unresolved; retain the unchanged KH clinical record, do not silently correct age.')
            rows.append({'patient_id':patient,'cohort':cohort['name'],'case_code':expected,'raw_mngs_patient_id':f'P{gid}','specimen_code':entry['specimen_id'],'raw_collection_date':raw_data[0]['collected_time'],'raw_mngs_sha256':sha(raw),'clinical_files_sha256':hashes,'notes':notes})
    return {'schema_version':'direct_raw_benchmark.input_mapping.v2','status':'frozen','source_index':{'file':'../data/patient_info_standardized/legacy_patient_index.json','sha256':sha(index_file)},'mapping_basis':'Reviewed case-code/specimen-level identity mapping. No Gold or model predictions used to generate this mapping.','patients':rows}

if __name__=='__main__':
    target=BASE/'input_mapping_v2.json'
    if target.exists(): raise SystemExit('Frozen mapping already exists; create a new version instead of overwriting it.')
    target.write_text(json.dumps(build_mapping(),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('Frozen 41 explicit input mappings: '+str(target))
