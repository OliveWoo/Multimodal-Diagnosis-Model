"""Replay recovered v19 against explicit, SHA256-pinned inputs without APIs.

This starts at saved ranked mNGS + deterministic summary, not raw/LLM generation.
The optional frozen merge is checked at its deterministic handoff; its advisory
LLM review and taxonomy annotations are preserved as frozen inputs, not rebuilt.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = 'research.upstream_v19_replay.v1'
PROFILE = 'docs/UPSTREAM_RECOVERY_PROFILE.json'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def object_json(data, label):
    value = json.loads(data.decode('utf-8-sig'))
    if not isinstance(value, dict):
        raise ValueError(f'{label}: expected a JSON object')
    return value


def verify_profile(root=ROOT):
    raw = (root/PROFILE).read_bytes()
    profile = object_json(raw, 'profile')
    if profile.get('profile_id') != 'recovered-v19-20260918':
        raise ValueError('Unexpected recovery profile')
    for row in profile['files']:
        path = (root/row['path']).resolve()
        if not path.is_relative_to(root.resolve()) or sha(path.read_bytes()) != row['sha256']:
            raise ValueError(f"Code/rules changed: {row['path']}")
    return sha(raw)


def read_reference(reference, base):
    if not isinstance(reference, dict) or not isinstance(reference.get('path'), str):
        raise ValueError('Every input requires an explicit path and SHA256')
    digest = reference.get('sha256', '')
    if not isinstance(digest, str) or not re.fullmatch(r'[0-9a-f]{64}', digest):
        raise ValueError('Missing or malformed input SHA256')
    path = (base/reference['path']).resolve()
    data = path.read_bytes()
    if sha(data) != digest:
        raise ValueError('Input SHA256 mismatch; refusing to select a replacement')
    return object_json(data, 'input'), path


def sans_source_paths(value):
    # Only these provenance containers are excluded. Rule versions, taxonomy,
    # array ordering, reasons, levels and all other content remain in comparison.
    if isinstance(value, dict):
        return {k: sans_source_paths(v) for k, v in value.items()
                if k not in {'source_files', 'merged_source_files'}}
    if isinstance(value, list):
        return [sans_source_paths(v) for v in value]
    return value


def deterministic_handoff(value):
    name = lambda x: str(x.get('organism_name', '')).strip().casefold()
    picked = [name(x) for x in value['best_available_summary']['picked_pathogens']]
    levels = {name(x): x.get('integrated_causative_level') for x in value['pathogen_candidates']}
    return {'picked_in_order': picked, 'candidate_levels': levels}


def ranked_for_patient(payload, patient_id, mngs):
    # Same payload selection/deduplication as mngs.load_ranked_mngs_for_patient,
    # operating on the already hash-checked bytes; no filename fallback.
    if isinstance(payload.get('records'), list):
        actual_id = str(payload.get('patient_id') or patient_id)
        selected = payload
    else:
        patients = payload.get('patients')
        if not isinstance(patients, dict):
            raise ValueError('Ranked input has neither records nor patients')
        selected = patients.get(patient_id) or patients.get(f'NGS_patient_{patient_id}')
        if not isinstance(selected, dict):
            raise ValueError('Requested patient absent from ranked input')
        actual_id = patient_id
    if actual_id.removeprefix('NGS_patient_') != patient_id:
        raise ValueError('Ranked patient ID mismatch')
    return mngs.deduplicate_ranked_mngs_payload({
        'patient_id': actual_id, 'records': selected.get('records', []),
        'ranking_metadata': selected.get('ranking_metadata', {}),
    })


def replay(manifest_path, output):
    manifest_path = manifest_path.resolve()
    output = output.resolve()
    if output.exists():
        raise ValueError('Output already exists; choose a new directory')
    profile_hash = verify_profile()
    manifest_raw = manifest_path.read_bytes()
    manifest = object_json(manifest_raw, 'manifest')
    if manifest.get('schema') != SCHEMA or not isinstance(manifest.get('cases'), list) or not manifest['cases']:
        raise ValueError('Expected a nonempty upstream_v19 replay manifest')
    prepared = []
    seen = set()
    for row in manifest['cases']:
        key = row.get('case_key', '')
        pid = row.get('patient_id', '')
        if not isinstance(key, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,80}', key) or key.casefold() in seen:
            raise ValueError('Invalid or duplicate case key')
        if not isinstance(pid, str) or not re.fullmatch(r'[0-9]+', pid):
            raise ValueError('Patient ID must be an explicit numeric string')
        seen.add(key.casefold())
        payloads = {}
        for label in ['ranked', 'summary', 'expected_scorer', 'frozen_merge']:
            if label == 'frozen_merge' and label not in row:
                continue
            payloads[label], path = read_reference(row.get(label), manifest_path.parent)
            if path.is_relative_to(output):
                raise ValueError('Output would contain an input file')
        prepared.append((key, pid, payloads, row))
    sys.path[:0] = [str(ROOT/'scripts/offline_guard'), str(ROOT/'upstream')]
    import sitecustomize  # block sockets and remove API keys before imports
    from tools import deterministic_mngs_max_scorer_v19_recovered as scorer
    results = []
    generated = []
    for key, pid, inputs, row in prepared:
        ranked = ranked_for_patient(inputs['ranked'], pid, scorer.mngs)
        fresh = scorer.score_payload(patient_dir=Path(f'NGS_patient_{pid}_json'),
                                     ranked_mngs=ranked, final_summary=inputs['summary'])
        result = {'case_key': key, 'scorer_content_equal':
                  sans_source_paths(fresh) == sans_source_paths(inputs['expected_scorer'])}
        if 'frozen_merge' in inputs:
            result['frozen_merge_handoff_equal'] = (
                deterministic_handoff(fresh) == deterministic_handoff(inputs['frozen_merge']['deterministic_max']))
        result['input_sha256'] = {k: row[k]['sha256'] for k in inputs}
        results.append(result)
        generated.append((key, fresh))
    receipt = {
        'profile_id': 'recovered-v19-20260918', 'profile_sha256': profile_hash,
        'manifest_sha256': sha(manifest_raw), 'cases': len(results),
        'scorer_content_equal': sum(r['scorer_content_equal'] for r in results),
        'frozen_merge_checked': sum('frozen_merge_handoff_equal' in r for r in results),
        'frozen_merge_handoff_equal': sum(r.get('frozen_merge_handoff_equal', False) for r in results),
        'comparison_exclusions': ['source_files', 'merged_source_files'],
        'raw_to_summary_rebuilt': False, 'llm_review_regenerated': False,
        'clinical_accuracy_evaluated': False, 'results': results,
    }
    # All inputs are checked and computation completed before creating output.
    output.mkdir(parents=True, exist_ok=False)
    for key, fresh in generated:
        with (output/f'{key}_scorer.json').open('x', encoding='utf-8') as f:
            json.dump(fresh, f, ensure_ascii=False, indent=2)
    with (output/'receipt.json').open('x', encoding='utf-8') as f:
        json.dump(receipt, f, ensure_ascii=False, indent=2)
    success = all(r['scorer_content_equal'] and r.get('frozen_merge_handoff_equal', True) for r in results)
    print(json.dumps({k: v for k, v in receipt.items() if k != 'results'}, ensure_ascii=True))
    return 0 if success else 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    try:
        return replay(args.manifest, args.output)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f'Replay refused: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
