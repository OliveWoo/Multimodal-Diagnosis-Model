"""Offline Excel/PPT frontend. Original research modules remain unchanged.

This adapter explicitly binds source files, frozen vision responses and patient
identity. It never calls a model. Run in its own process to isolate `tools`.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import logging
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SHARED = ROOT / 'reference_versions/shared_upstream'


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def semantic_sha(value):
    return hashlib.sha256(json.dumps(value, ensure_ascii=False, sort_keys=True,
                                     separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def safe_name(value):
    require(isinstance(value, str) and bool(re.fullmatch(r'[A-Za-z0-9_.-]+', value))
            and value not in {'.', '..'}, 'Unsafe output name')
    return value


def identity_rows(rows):
    return [[r['case_code'], r['specimen_id'], r['hospital_id']] for r in rows]


def validate_identity(rows, expected=None):
    identities = identity_rows(rows)
    require(bool(rows), 'Empty patient index')
    require(len({tuple(r) for r in identities}) == len(rows), 'Duplicate specimen identity')
    require(len({r['patient_key'] for r in rows}) == len(rows), 'Duplicate patient key')
    for r in rows:
        require(r['patient_key'] == Path(r['patient_key']).name and '/' not in r['patient_key']
                and '\\' not in r['patient_key'] and r['patient_key'] not in {'.', '..'}, 'Unsafe patient key')
    digest = semantic_sha(identities)
    if expected:
        require(digest == expected, 'Patient identity/order differs from locked index')
    return digest


def initialize():
    sys.dont_write_bytecode = True
    # Explicit execution also works when Python started without PYTHONPATH.
    import runpy
    runpy.run_path(str(ROOT / 'scripts/offline_guard/sitecustomize.py'))
    require('tools' not in sys.modules, 'Start frontend in a separate Python process')
    sys.path.insert(0, str(SHARED))


class Inputs:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.used = {}

    def file(self, ref):
        require(isinstance(ref, dict) and set(ref) == {'path', 'sha256'}, 'Expected path + sha256 reference')
        path = (self.directory / ref['path']).resolve()
        require(path.is_file(), 'Referenced input file is missing')
        require(sha(path) == ref['sha256'], 'Input SHA256 mismatch: ' + path.name)
        self.used[str(path)] = ref['sha256']
        return path

    def unchanged(self):
        require(all(Path(p).is_file() and sha(p) == h for p, h in self.used.items()), 'An input changed during execution')


def verify_profile(expected):
    profile_path = ROOT / 'docs/FRONTEND_CODE_PROFILE.json'
    require(sha(profile_path) == expected, 'Frontend code profile differs from manifest')
    for row in read(profile_path)['files']:
        require(sha(ROOT / row['path']) == row['sha256'], 'Frontend module changed: ' + row['path'])


def strip_image_paths(value):
    if isinstance(value, dict):
        return {k: strip_image_paths(v) for k, v in value.items() if k != 'raw_path'}
    if isinstance(value, list):
        return [strip_image_paths(v) for v in value]
    return value


class FrozenVision:
    model = 'frozen-response-replay'
    fallback_models = ()

    def __init__(self, source, bindings, inputs):
        self.expected = {}
        self.records = {}
        self.used = set()
        for ref in bindings:
            path = inputs.file(ref)
            payload = read(path)
            require(payload['source_file'] == source.name, 'Frozen slide source filename mismatch')
            slide = payload['slide_index']
            require(isinstance(slide, int) and slide > 0 and slide not in self.expected, 'Duplicate/invalid frozen slide')
            self.expected[slide] = payload
            for record in payload['images']:
                key = (slide, record['image_index'])
                require(key not in self.records, 'Duplicate frozen image binding')
                filename = Path(record['raw_path'].replace('\\', '/')).name
                image = path.parent / 'images' / filename
                # Frozen parsed.json is hash-pinned; its sibling image bytes are
                # checked against the image extracted from the pinned PPT.
                require(image.is_file(), 'Frozen image file missing')
                digest = sha(image)
                inputs.used[str(image.resolve())] = digest
                vision = record['vision']
                require(vision.get('status') in {'success', 'empty'}, 'Unsuccessful frozen image response')
                require(all(isinstance(vision.get(k), list) and all(isinstance(s, str) for s in vision[k])
                            for k in ('raw_text', 'findings')), 'Invalid frozen response text')
                self.records[key] = (digest, vision)

    def extract(self, image_path):
        from ppt_patient_parser.vision_hooks import VisionExtractionResult
        key = (int(image_path.parent.parent.name.split('_')[-1]), int(image_path.stem.split('_')[-1]))
        require(key in self.records, 'No frozen response for this slide/image; network is disabled')
        digest, value = self.records[key]
        require(sha(image_path) == digest, 'PPT image bytes differ from frozen response image')
        require(key not in self.used, 'Image response reused within a slide')
        self.used.add(key)
        return VisionExtractionResult(raw_text=value['raw_text'], findings=value['findings'],
                                      confidence=value.get('confidence'), model_used=value.get('model_used'))


def extract_ppt(item, inputs, target):
    from pptx import Presentation
    from ppt_patient_parser.cli import process_presentation
    source = inputs.file(item['source'])
    frozen = FrozenVision(source, item.get('frozen_slides', []), inputs)
    process_presentation(source, target, logging.getLogger('frontend'), vision_extractor=frozen)
    parsed = sorted(target.rglob('parsed.json'))
    require(len(parsed) == len(Presentation(source).slides), 'PPT parser skipped one or more slides')
    require(bool(parsed), 'PPT has no slides')
    comparisons = 0
    for path in parsed:
        value = read(path)
        require(not read(path.with_name('debug.json'))['errors'], 'PPT extraction reported an error; see private log')
        if frozen.expected:
            require(value['slide_index'] in frozen.expected, 'Frozen slide is missing')
            require(strip_image_paths(value) == strip_image_paths(frozen.expected[value['slide_index']]),
                    'Parsed slide differs from saved result')
            comparisons += 1
    require(frozen.used == set(frozen.records), 'Unused frozen image bindings')
    require(not frozen.expected or len(frozen.expected) == len(parsed), 'Extra frozen slides')
    return parsed, {'slides': len(parsed), 'images_replayed': len(frozen.used), 'saved_slides_equal': comparisons}


@contextlib.contextmanager
def radiology_route(linker, mode):
    require(mode in {'labs_only', 'radiology_and_labs'}, 'Unknown PPT import route')
    original = linker.extract_radiology_entries
    try:
        if mode == 'labs_only':
            linker.extract_radiology_entries = lambda payload: []
        yield
    finally:
        linker.extract_radiology_entries = original


def master_ppt(config, inputs, output, receipt):
    from tools.split_excel_master_by_patient import export_patient_workbook
    from tools.export_legacy_patient_json import export_legacy_patient_json
    from tools import link_ppt_slides_to_excel_patients as linker
    from tools.import_ppt_labs_to_patient_exports import import_ppt_labs
    import pandas as pd
    source = inputs.file(config['master'])
    with pd.ExcelFile(source) as book:
        require({'基本資料', '微生物診斷資料', 'mNGS', 'Read count'} <= set(book.sheet_names), 'Unsupported master workbook schema')
        for sheet, columns in {
            '基本資料': {'編碼', 'Specimen_ID', 'Hospital_ID'},
            '微生物診斷資料': {'編碼', 'Specimen_ID', 'Hospital_ID'},
            'mNGS': {'編碼', 'specimen_id', 'hospital_patient_id'},
            'Read count': {'specimen_id', 'hospital_patient_id'},
        }.items():
            require(columns <= set(pd.read_excel(book, sheet_name=sheet, nrows=0).columns),
                    'Required identity columns missing in ' + sheet)
    patients = output / 'patient_exports'
    _, index = export_patient_workbook(source, patients)
    rows = read(index)
    receipt['identity_sha256'] = validate_identity(rows, config.get('expected_identity_sha256'))
    for row in rows:
        require(Path(row['output_dir']).resolve() == (patients / row['patient_key']).resolve(), 'Patient output path escapes run')
    patient_index, by_hospital = linker._load_patient_index(patients)
    digests, names, bindings, batches = set(), set(), [], []
    for i, item in enumerate(config.get('presentations', []), 1):
        ppt = inputs.file(item['source'])
        require(item['source']['sha256'] not in digests, 'Duplicate PPT content; deduplicate explicitly in manifest')
        require(ppt.name not in names, 'PPT filename collision would merge provenance')
        digests.add(item['source']['sha256']); names.add(ppt.name)
        target = output / 'ppt' / f'input_{i:03}'
        parsed, counts = extract_ppt(item, inputs, target)
        for path in parsed:
            payload = read(path)
            keys, reason = linker._resolve_patients_for_slide(payload, patient_index, by_hospital)
            bindings.append({'ppt_sha256': item['source']['sha256'], 'slide_index': payload['slide_index'],
                             'patient_keys': keys, 'reason': reason, 'route': item['route']})
        with radiology_route(linker, item['route']):
            linked = linker.link_ppt_to_patient_exports(target, patients)
        batches.append(dict(source_sha256=item['source']['sha256'], route=item['route'], **counts, link=linked))
    binding_digest = semantic_sha(bindings)
    if config.get('expected_binding_sha256'):
        require(binding_digest == config['expected_binding_sha256'], 'PPT-to-patient binding changed')
    write(output / 'ppt_bindings.private.json', bindings)
    write(output / 'ppt_batches.private.json', batches)
    receipt.update(cases=len(rows), presentations=len(batches), slides=sum(b['slides'] for b in batches),
                   frozen_images=sum(b['images_replayed'] for b in batches), binding_sha256=binding_digest,
                   unmatched_slides=sum(not b['patient_keys'] for b in bindings),
                   fanout_slides=sum(len(b['patient_keys']) > 1 for b in bindings))
    receipt['labs'] = import_ppt_labs(patients)
    export_legacy_patient_json(patients, output / 'standardized')
    receipt['stage'] = 'standardized_with_frozen_vision' if receipt['frozen_images'] else 'standardized'


def legacy_excel(config, inputs, output, receipt):
    from core.raw_parser import convert_workbook
    from core.llm_parser import _build_prompt
    source = inputs.file(config['workbook'])
    succeeded, errors, _ = convert_workbook(source, output_dir=output / 'raw')
    require(not errors and bool(succeeded), 'Excel sheet extraction failed; see private log')
    raw_files = sorted((output / 'raw').glob('*.json'))
    # The shared source contains exact prompt names for all eight raw sheets.
    # Resolve against the pinned source tree, never the caller's current folder.
    requests = []
    for path in raw_files:
        payload = read(path)
        prompt = SHARED / 'prompts' / (payload['sheet'] + '_prompt.txt')
        require(prompt.resolve().is_relative_to((SHARED / 'prompts').resolve()) and prompt.is_file(),
                'No exact formatter prompt for raw sheet')
        requests.append({'raw_name': path.name, 'raw_semantic_sha256': semantic_sha(payload),
                         'prompt_file': prompt.relative_to(ROOT).as_posix(), 'prompt_sha256': sha(prompt),
                         'historical_shared_code_model': 'gpt-5', 'historical_shared_code_temperature': 1,
                         'request_text': _build_prompt(prompt.read_text(encoding='utf-8'), payload)})
    write(output / 'formatter_requests.private.json', requests)
    receipt.update(raw_sections=len(raw_files), prepared_formatter_requests=len(requests),
                   stage='raw_extracted_formatter_prepared_not_run')
    responses = config.get('frozen_formatter')
    if responses is not None:
        require(bool(config.get('identity_namespace')), 'Frozen formatter replay requires an explicit patient namespace')
        covered, destinations = set(), set()
        for item in responses:
            name = safe_name(item['raw_name'])
            destination = safe_name(item['output_name'])
            require(name not in covered and destination not in destinations, 'Duplicate formatter binding')
            require(destination.endswith('.json'), 'Formatter output must be JSON')
            covered.add(name); destinations.add(destination)
            require(semantic_sha(read(output / 'raw' / name)) == item['raw_semantic_sha256'], 'Raw data differs from frozen formatter input')
            response = inputs.file(item['response'])
            require(isinstance(read(response), (list, dict)), 'Frozen formatter output is not structured JSON')
            write(output / 'standardized' / destination, read(response))
            if item.get('prompt'):
                inputs.file(item['prompt'])
        require(covered == {p.name for p in raw_files}, 'Incomplete frozen formatter coverage')
        receipt.update(stage='normalized_from_explicit_frozen_responses', frozen_formatter_responses=len(responses),
                       historical_prompt_binding='User-supplied bindings; not proof of historical model execution')


def table1(config, inputs, output, receipt):
    from tools.parse_table1_standardized import export_standardized, STANDARD_SUFFIXES
    from tools.parse_table1_workbook import parse_workbook
    source = inputs.file(config['workbook'])
    patient_id = safe_name(config['patient_id'])
    require(bool(config.get('identity_namespace')), 'Table 1 requires an explicit patient namespace')
    base = patient_id if patient_id.startswith('NGS_patient_') else 'NGS_patient_' + patient_id
    baseline = None
    if config.get('baseline') is not None:
        require(bool(config.get('baseline_identity_evidence')), 'Explain baseline identity; numeric patient ID alone is insufficient')
        baseline = output / 'baseline_inputs'
        for suffix in STANDARD_SUFFIXES:
            path = inputs.file(config['baseline'][suffix])
            baseline.mkdir(exist_ok=True)
            shutil.copyfile(path, baseline / f'{base}_{suffix}.json')
    sheet = config.get('sheet', 'Table 1')
    parse_workbook(source, output / 'raw', sheet_name=sheet)
    directory = export_standardized(source, output / 'standardized', patient_id=patient_id,
                                    sheet_name=sheet, baseline_dir=baseline)
    receipt.update(stage='standardized', clinical_sections=len(STANDARD_SUFFIXES),
                   identity_namespace=config['identity_namespace'], baseline_used=baseline is not None,
                   standardized_directory=str(directory.relative_to(output)))


def compare_expected(expected, inputs, output):
    checks = []
    for row in expected:
        actual = (output / row['output']).resolve()
        require(actual.is_relative_to(output.resolve()), 'Expected output escapes run directory')
        reference = inputs.file(row['reference'])
        equal = actual.is_file() and read(actual) == read(reference)
        checks.append({'output': row['output'], 'equal': equal})
    write(output / 'comparison.private.json', checks)
    require(all(r['equal'] for r in checks), 'Output content differs from expected files')
    return len(checks)


def run(manifest, output):
    manifest = Path(manifest).resolve(); output = Path(output).resolve()
    require(not output.exists(), 'Output already exists; choose a new directory')
    config = read(manifest)
    require(config.get('schema_version') == 1, 'Unsupported manifest schema')
    mode = config.get('mode')
    require(mode in {'master_ppt', 'legacy_excel', 'table1'}, 'Unsupported frontend mode')
    require(config.get('run_kind') in {'extract', 'replay'}, 'Choose extract or replay explicitly')
    if config['run_kind'] == 'replay':
        require(bool(config.get('expected')), 'Replay requires explicit saved-output comparisons')
        if mode == 'master_ppt':
            require(bool(config.get('expected_identity_sha256')) and bool(config.get('expected_binding_sha256')),
                    'Master replay requires patient identity/order and PPT binding locks')
    verify_profile(config['code_profile_sha256'])
    # Validate all path references before creating any output.
    inputs = Inputs(manifest.parent)
    def visit(value):
        if isinstance(value, dict):
            if set(value) == {'path', 'sha256'}:
                inputs.file(value)
            else:
                for child in value.values(): visit(child)
        elif isinstance(value, list):
            for child in value: visit(child)
    visit(config)
    require(all(not Path(p).is_relative_to(output) for p in inputs.used), 'Output contains a source input')
    output.mkdir(parents=True)
    receipt = dict(schema_version=1, mode=mode, run_kind=config['run_kind'], status='running',
                   manifest_sha256=sha(manifest), code_profile_sha256=config['code_profile_sha256'],
                   network='disabled', fresh_model_calls=0)
    logging.basicConfig(filename=output / 'frontend.private.log', level=logging.ERROR, encoding='utf-8', force=True)
    try:
        globals()[mode](config, inputs, output, receipt)
        receipt['exact_json_comparisons'] = compare_expected(config.get('expected', []), inputs, output)
        inputs.unchanged()
        receipt['status'] = 'passed'
    except Exception as error:
        receipt.update(status='failed', error=str(error))
        raise
    finally:
        receipt['inputs'] = [{'path': p, 'sha256': h} for p, h in sorted(inputs.used.items())]
        write(output / 'receipt.private.json', receipt)
    return {k: v for k, v in receipt.items() if k not in {'inputs', 'labs'}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    initialize()
    try:
        print(json.dumps(run(args.manifest, args.output), ensure_ascii=True))
    except Exception as error:
        print(json.dumps({'status': 'failed', 'error': str(error)}, ensure_ascii=True), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
