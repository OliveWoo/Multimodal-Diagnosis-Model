"""Keep KH case-code IDs and two-hospital global IDs in separate handoffs.

This entry validates and stages explicitly frozen clinical inputs. It is not a
replacement for Excel parsing, model inference, or mNGS candidate generation.
"""
import argparse
import json
import re
import shutil
import sys
from pathlib import Path

from pipeline import ROOT, Inputs, initialize, read, require, sha, write

CLINICAL = ('admission_diagnosis', 'underlying', 'CBC', 'other_lab', 'culture',
            'filmarray', 'gm_test', 'image')
PROFILES = {
    'KH': {'namespace': 'kmuh_case_code', 'sections': CLINICAL,
           'anchor': 'all_RK_NTC_microbes', 'prefixes': {'K'}},
    'two_hospital': {'namespace': 'legacy_global', 'sections': CLINICAL + ('molecular_microbiology',),
                     'anchor': 'mNGS_grouped', 'prefixes': {'T', 'V'}},
}


def normalized(value):
    return re.sub(r'[^A-Z0-9]', '', str(value or '').upper())


def patient_number(value):
    require(isinstance(value, str) and bool(re.fullmatch(r'P[1-9]\d*', value)), 'Use an explicit P<number> identifier')
    return value[1:]


def validate_binding(cohort, case, index):
    local_id = patient_number(case['patient_id'])
    global_id = case['master_patient_id']
    patient_number(global_id)
    require(global_id in index, 'Master patient ID is missing from pinned index')
    row = index[global_id]
    require(case['case_code'] == row['case_code'], 'Case code differs from pinned index')
    require(normalized(case['specimen_code']) == normalized(row['specimen_id']), 'Specimen differs from pinned index')
    require(str(case['case_code'])[:1] in PROFILES[cohort]['prefixes'], 'Hospital group is incompatible with this entry')
    if cohort == 'KH':
        require(case['case_code'] == f'K{int(local_id):03}', 'KH clinical ID must follow K case code, not global ordinal')
    else:
        require(case['patient_id'] == global_id, 'Two-hospital global IDs must not be renumbered')
    return row


def validate_anchor(cohort, case, value):
    require(isinstance(value, list), 'mNGS anchor must be a specimen-record list')
    if not value:
        require(cohort == 'two_hospital' and case.get('anchor_status') == 'empty_explicit',
                'Empty specimen anchor must be explicitly recorded; KH requires its specimen record')
        return 'empty_mngs_anchor_identity_from_index_only'
    require(case.get('anchor_status') == 'specimen_checked', 'Declare specimen_checked for a populated anchor')
    require(all(isinstance(r, dict) and isinstance(r.get('pathogens'), dict) and r.get('specimen_code') for r in value),
            'Flat RK/NTC candidate lists are not grouped specimen records')
    require({normalized(r['specimen_code']) for r in value} == {normalized(case['specimen_code'])},
            'mNGS specimen does not match this clinical case')
    return 'specimen_checked'


def check_filename(path, patient_id, section):
    # A numbered OS duplicate suffix is accepted only after explicit hash binding.
    clean = re.sub(r'\s+\(\d+\)(?=\.json$)', '', path.name)
    require(clean == f'NGS_patient_{patient_number(patient_id)}_{section}.json',
            'Clinical filename belongs to a different local patient/section')


def verify_code(expected):
    profile = ROOT / 'docs/COHORT_CODE_PROFILE.json'
    require(sha(profile) == expected, 'Cohort code profile mismatch')
    for item in read(profile)['files']:
        require(sha(ROOT / item['path']) == item['sha256'], 'Cohort code changed: ' + item['path'])


def run(cohort, manifest, output):
    require(cohort in PROFILES, 'Unsupported cohort')
    manifest = Path(manifest).resolve(); output = Path(output).resolve()
    require(not output.exists(), 'Output already exists; refusing overwrite')
    config = read(manifest)
    require(config.get('schema_version') == 'research_frontend.cohort.v1', 'Unsupported cohort manifest')
    require(config.get('cohort') == cohort, 'Manifest is for the other hospital group')
    require(config.get('clinical_id_namespace') == PROFILES[cohort]['namespace'], 'Wrong clinical ID namespace')
    require(config.get('operation') == 'prepare_verified_snapshot', 'This entry stages saved inputs only')
    require(config.get('source_boundary') == 'saved_clinical_and_mngs', 'Declare the frozen source boundary')
    require(bool(config.get('dataset_id')), 'Declare a source dataset version')
    verify_code(config['code_profile_sha256'])
    inputs = Inputs(manifest.parent)
    def visit(v):
        if isinstance(v, dict):
            if set(v) == {'path', 'sha256'}: inputs.file(v)
            else:
                for item in v.values(): visit(item)
        elif isinstance(v, list):
            for item in v: visit(item)
    visit(config)
    index = {}
    for row in read(inputs.file(config['source_index'])):
        global_id = row['legacy_patient_id'].replace('NGS_patient_', 'P')
        patient_number(global_id)
        require(global_id not in index, 'Duplicate global ID in source index')
        index[global_id] = row
    cases = config['cases']
    require(bool(cases), 'Empty cohort')
    local_ids, global_ids, prepared = set(), set(), []
    for case in cases:
        require(case['patient_id'] not in local_ids, 'Duplicate clinical patient ID')
        require(case['master_patient_id'] not in global_ids, 'Duplicate master identity in cohort')
        local_ids.add(case['patient_id']); global_ids.add(case['master_patient_id'])
        validate_binding(cohort, case, index)
        require(set(case['clinical']) == set(PROFILES[cohort]['sections']), 'Clinical section coverage differs from cohort contract')
        anchor_name = PROFILES[cohort]['anchor']
        anchor = inputs.file(case['mngs_anchor'])
        check_filename(anchor, case['patient_id'], anchor_name)
        identity_status = validate_anchor(cohort, case, read(anchor))
        sources = dict(case['clinical'], **{anchor_name: case['mngs_anchor']})
        expected = case.get('expected_upstream', {})
        require(not expected or set(expected) == set(sources), 'Upstream comparison must cover all staged sections')
        for section, ref in sources.items():
            path = inputs.file(ref); check_filename(path, case['patient_id'], section)
            require(isinstance(read(path), (list, dict)), 'Clinical input is not structured JSON')
        prepared.append((case, sources, expected, identity_status))
    for name, selected in config.get('subsets', {}).items():
        require(bool(re.fullmatch(r'[A-Za-z0-9_-]+', name)), 'Unsafe subset name')
        require(len(selected) == len(set(selected)) and set(selected) <= local_ids, 'Subset contains duplicate/unknown patients')
    require(all(not Path(p).is_relative_to(output) for p in inputs.used), 'Output contains an input file')
    output.mkdir(parents=True)
    cohort_root = output / cohort
    receipt = dict(cohort=cohort, dataset_id=config['dataset_id'], clinical_id_namespace=config['clinical_id_namespace'],
                   operation=config['operation'], source_boundary=config['source_boundary'], status='running',
                   manifest_sha256=sha(manifest), code_profile_sha256=config['code_profile_sha256'],
                   fresh_model_calls=0, full_raw_to_ober_reproduced=False,
                   patient_count=len(cases), subset_counts={k:len(v) for k,v in config.get('subsets', {}).items()})
    mappings, comparisons, files = [], [], []
    try:
        for case, sources, expected, identity_status in prepared:
            base = 'NGS_patient_' + patient_number(case['patient_id'])
            target = cohort_root / 'patient_info' / (base + '_json')
            target.mkdir(parents=True)
            for section, ref in sources.items():
                source = inputs.file(ref); destination = target / f'{base}_{section}.json'
                shutil.copyfile(source, destination)
                require(sha(source) == sha(destination), 'Staged input bytes changed')
                files.append({'patient_id':case['patient_id'], 'section':section,
                              'output':destination.relative_to(output).as_posix(), 'sha256':sha(destination)})
                if expected:
                    comparisons.append({'patient_id':case['patient_id'], 'section':section,
                                        'equal':read(destination) == read(inputs.file(expected[section]))})
            mappings.append({k:case[k] for k in ('patient_id','master_patient_id','case_code','specimen_code')} |
                            {'clinical_id_namespace':config['clinical_id_namespace'], 'identity_status':identity_status})
        write(cohort_root / 'identity_mapping.private.json', mappings)
        write(cohort_root / 'files.private.json', files)
        write(cohort_root / 'upstream_comparison.private.json', comparisons)
        write(cohort_root / 'subsets.private.json', config.get('subsets', {}))
        receipt.update(staged_files=len(files), upstream_sections_compared=len(comparisons),
                       upstream_sections_equal=sum(r['equal'] for r in comparisons),
                       empty_mngs_anchor_cases=sum(r['identity_status'].startswith('empty_') for r in mappings))
        require(all(r['equal'] for r in comparisons), 'Staged data differs from the selected upstream snapshot')
        inputs.unchanged()
        receipt['status'] = 'passed'
    except Exception as error:
        receipt.update(status='failed', error=str(error))
        raise
    finally:
        receipt['inputs'] = [{'path':p,'sha256':h} for p,h in sorted(inputs.used.items())]
        write(output / 'receipt.private.json', receipt)
    return {k:v for k,v in receipt.items() if k != 'inputs'}


def main(cohort):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    initialize()
    try:
        print(json.dumps(run(cohort, args.manifest, args.output), ensure_ascii=True))
        return 0
    except Exception as error:
        print(json.dumps({'status':'failed','error':str(error)}, ensure_ascii=True), file=sys.stderr)
        return 1
