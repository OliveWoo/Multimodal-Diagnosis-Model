"""Local-only source-identity validation, deliberately separate from model payloads."""
import hashlib
import json
import re
from pathlib import Path

def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def read_json(path): return json.loads(path.read_text(encoding='utf-8-sig'))
def normalized_code(v): return re.sub(r'[^A-Z0-9]','',str(v or '').upper())

def read_frozen_json(path,expected_hash):
    data=path.read_bytes()
    if hashlib.sha256(data).hexdigest()!=expected_hash:
        raise ValueError(f'Source hash mismatch: {path.name}')
    return json.loads(data.decode('utf-8-sig'))

def load_mapping(manifest,base,sections):
    if not manifest.get('input_mapping'):
        raise ValueError('Explicit input_mapping is required; implicit same-number joins are forbidden.')
    mapping_path=(base/manifest['input_mapping']).resolve()
    mapping=read_json(mapping_path)
    if mapping.get('schema_version')!='direct_raw_benchmark.input_mapping.v2' or mapping.get('status')!='frozen':
        raise ValueError('Require a frozen v2 input mapping, not an audit proposal.')
    source=mapping['source_index']
    source_path=(mapping_path.parent/source['file']).resolve()
    index=read_frozen_json(source_path,source['sha256'])
    index_by_id={}
    for row in index:
        global_id=row['legacy_patient_id'].replace('NGS_patient_','P')
        if global_id in index_by_id: raise ValueError('Duplicate global ID in source index')
        index_by_id[global_id]=row
    expected={(c['name'],p.upper()) for c in manifest['cohorts'] for p in c['patient_ids']}
    namespaces={c['name']:c.get('clinical_id_namespace') for c in manifest['cohorts']}
    result={}
    for row in mapping['patients']:
        key=(row['cohort'],row['patient_id'])
        if key in result: raise ValueError('Duplicate mapping row')
        if key not in expected: raise ValueError('Mapping contains a patient outside the frozen cohort')
        if not re.fullmatch(r'P[1-9]\d*',row['raw_mngs_patient_id']): raise ValueError('Invalid raw mNGS patient ID')
        source_row=index_by_id.get(row['raw_mngs_patient_id'])
        if not source_row or source_row['case_code']!=row['case_code'] or normalized_code(source_row['specimen_id'])!=normalized_code(row['specimen_code']):
            raise ValueError('Mapping case/specimen does not match the source index')
        namespace=namespaces[row['cohort']]
        if namespace=='kmuh_case_code':
            if row['case_code']!=f"K{int(row['patient_id'][1:]):03}": raise ValueError('KH clinical case code mismatch')
        elif namespace=='legacy_global':
            if row['raw_mngs_patient_id']!=row['patient_id'] or row['case_code'].startswith('K'):
                raise ValueError('Two-hospital legacy identity mismatch')
        else: raise ValueError('Declare the clinical ID namespace explicitly')
        if set(row['clinical_files_sha256'])!=set(sections): raise ValueError('Clinical source hash coverage is incomplete')
        for h in [row['raw_mngs_sha256'],*row['clinical_files_sha256'].values()]:
            if h is not None and not re.fullmatch(r'[0-9a-f]{64}',h): raise ValueError('Invalid SHA256 in mapping')
        result[key]=row
    if set(result)!=expected: raise ValueError('Input mapping does not cover the entire frozen cohort')
    return result,sha256_file(mapping_path)
