"""Offline inspection of the supplied v7 result release; never runs inference.

Only writes an isolated source copy and descriptive audit JSON. Evaluation labels
are read solely for this retrospective review, never for a proposed OBER input.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import stat
import zipfile
from collections import Counter, defaultdict
from pathlib import Path, PurePosixPath


HERE = Path(__file__).resolve().parent
ARCHIVE = Path('C:/Users/chung/Downloads/2026-10-20_integrated_v7_latest_results_no_OBER.zip')
SOURCE = HERE / 'source'
COHORTS = ['KH33', 'KH18', 'A11', 'B11']
POSITIVE_TIERS = {'Picked', 'Possible', 'Fallback-Possible'}

# Retrospective name reconciliation only, not a new production scoring contract.
# Complex/group mappings below are observed in this release's reference flags;
# they must NOT be reused as organism-specific clinical evidence in an OBER judge.
DIAGNOSTIC_NAME_MAP = {
    'cmv': 'humanbetaherpesvirus5',
    'humancytomegalovirus': 'humanbetaherpesvirus5',
    'hsv': 'humanalphaherpesvirus1',
    'hsv1': 'humanalphaherpesvirus1',
    'herpessimplexvirustype1': 'humanalphaherpesvirus1',
    'vzv': 'humanalphaherpesvirus3',
    'covid19': 'sarscov2',
    'severeacuterespiratorysyndromerelatedcoronavirus': 'sarscov2',
    'crkp': 'klebsiellapneumoniae',
    'klebsiellapneumoniaegroup': 'klebsiellapneumoniae',
    'acinetobactercalcoaceticusbaumanniicomplex': 'acinetobacterbaumannii',
    'candidaglabrata': 'nakaseomycesglabratus',
    'candidaalbican': 'candidaalbicans',
    'influenxaa': 'influenzaavirus',
    'legionellapneumophilasubsppneumophila': 'legionellapneumophila',
    'aspergillosis': 'aspergillusspp',
    'ecoli': 'escherichiacoli',
    'mrsa': 'staphylococcusaureus',
}


def diagnostic_name(value):
    key = re.sub('[^a-z0-9]', '', value.lower())
    return DIAGNOSTIC_NAME_MAP.get(key, key)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def flag(value):
    require(value in ('True', 'False'), f'Unexpected boolean: {value!r}')
    return value == 'True'


def keys(row):
    return row['cohort'], str(row['patient_id']), row['organism_name']


def write_new(path, data):
    if path.exists():
        require(path.is_file() and path.read_bytes() == data, f'Refuse different existing output: {path}')
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        stream.write(data)


def main():
    archive_hash = digest(ARCHIVE.read_bytes())
    with zipfile.ZipFile(ARCHIVE) as archive:
        require(archive.testzip() is None, 'ZIP CRC failure')
        infos = [info for info in archive.infolist() if not info.is_dir()]
        require(len({i.filename.casefold() for i in infos}) == len(infos), 'Duplicate archive paths')
        tops = {PurePosixPath(i.filename).parts[0] for i in infos}
        require(len(tops) == 1, 'Expected one release root')
        root = SOURCE / next(iter(tops))
        plan = []
        for info in infos:
            pp = PurePosixPath(info.filename)
            require(not pp.is_absolute() and '..' not in pp.parts and '\\' not in info.filename
                    and ':' not in info.filename and not stat.S_ISLNK(info.external_attr >> 16), 'Unsafe archive path')
            dst = (SOURCE / str(pp)).resolve()
            require(dst.is_relative_to(SOURCE.resolve()), 'Destination outside source copy')
            require(pp.suffix.lower() in {'.json', '.md', '.csv'}, 'Unexpected archive type')
            content = archive.read(info)
            require(not dst.exists() or dst.read_bytes() == content, f'Conflicting source copy: {dst}')
            plan.append((dst, content))
        for dst, content in plan:
            write_new(dst, content)
    manifest = load(root / 'package_manifest.json')
    for entry in manifest['files']:
        p = (root / entry['package_path']).resolve()
        require(p.is_relative_to(root.resolve()) and p.is_file(), 'Manifest destination invalid')
        require(digest(p.read_bytes()) == entry['source_sha256'] == entry['package_sha256'],
                f'Manifest hash mismatch: {entry["package_path"]}')

    evaluation = load(root / 'postfreeze_evaluation/summary.json')
    candidate_labels = rows(root / 'postfreeze_evaluation/all_candidate_before_after.csv')
    patient_labels = rows(root / 'postfreeze_evaluation/patient_before_after.csv')
    label_by_key = {keys(r): r for r in candidate_labels}
    require(len(label_by_key) == len(candidate_labels) == 570, 'Candidate label keys not unique')
    require(len({(r['cohort'], r['patient_id']) for r in patient_labels}) == len(patient_labels) == 73,
            'Patient label keys not unique')
    audit = {
        'archive': str(ARCHIVE), 'archive_sha256': archive_hash,
        'source_root': str(root), 'source_file_count': len(infos),
        'manifest_entries_verified': len(manifest['files']),
        'ober_included': manifest['ober_included'], 'api_calls': 0,
        'metrics_basis': 'Supplied postfreeze summary, reconciled to candidate flags; no synonym changes.',
        'cohorts': {}, 'tier_totals': {}, 'missed_in_context': [], 'selected_false_positives': [],
        'candidate_source_counts': {}, 'decision_evidence_modes': {},
        'schema_observations': {}, 'patients': [],
        'diagnostic_name_reconciliation': {
            'mapping': DIAGNOSTIC_NAME_MAP,
            'scope': 'Retrospective diagnostic check, not frozen production aliases or a new official score.',
            'official_alias_contents_in_package': False,
        },
    }
    observed_keys = set()
    tier_total = defaultdict(Counter)
    candidate_sources = Counter()
    evidence_modes = Counter()
    shape_keys = Counter()
    old_merge_keys = Counter()
    viral_counts = Counter()
    all_blockers = Counter()
    all_case_ids = defaultdict(list)
    for cohort in COHORTS:
        by_tier = defaultdict(Counter)
        total_candidates = total_selected = selected_tp = context_tp = 0
        distinct_in_pool_fn = distinct_outside_pool_fn = 0
        reconstructed_tp = reconstructed_fp = reconstructed_fn = 0
        profiles = Counter()
        cases = []
        for path in sorted((root / cohort / 'patient_results').glob('*.json')):
            payload = load(path)
            pid = str(payload['patient_id'])
            all_case_ids[pid].append(cohort)
            candidates = payload['all_candidates']
            packet_path = root / cohort / 'reporting/patient_outputs' / f'NGS_patient_{pid}_test_aware_possible_pathogen_shadow.json'
            packet = load(packet_path)
            require(digest(packet_path.read_bytes()) == payload['source_reporting_sha256'], 'Packet hash mismatch')
            require(packet['candidate_count'] == payload['candidate_count'] == len(candidates), 'Candidate count mismatch')
            packet_by_name = {c['organism_name']: c for c in packet['all_forwarded_candidates']}
            require(len(packet_by_name) == len(candidates), 'Packet candidate keys mismatch')
            profile = payload['decision_evidence_profile']['decision_evidence_mode']
            profiles[profile] += 1
            evidence_modes[profile] += 1
            shape_keys.update(payload.keys())
            for old_key in ['deterministic_max', 'llm_missed_candidate_review']:
                old_merge_keys[old_key] += int(old_key in payload or old_key in packet)
            case_selected = []
            case_tiers = Counter()
            for candidate in candidates:
                name = candidate['organism_name']
                key = (cohort, pid, name)
                require(key not in observed_keys, f'Duplicate candidate: {key}')
                observed_keys.add(key)
                label = label_by_key[key]
                selected = candidate['reporting']['selected_for_complete_report']
                tier = candidate['reporting']['tier']
                correct = flag(label['is_reference_answer'])
                require(type(selected) is bool and selected == (tier in POSITIVE_TIERS), 'Selection/tier mismatch')
                require(selected == flag(label['after_selected']) and tier == label['after_tier'], 'Evaluation/output mismatch')
                source_candidate = packet_by_name[name]
                if selected:
                    case_selected.append(name)
                case_tiers[tier] += 1
                by_tier[tier]['candidate_count'] += 1
                by_tier[tier]['reference_positive'] += int(correct)
                by_tier[tier]['reference_nonmatching'] += int(not correct)
                tier_total[tier].update({'candidate_count': 1, 'reference_positive': int(correct),
                                         'reference_nonmatching': int(not correct)})
                total_candidates += 1
                total_selected += int(selected)
                selected_tp += int(selected and correct)
                context_tp += int(not selected and correct)
                candidate_sources[candidate.get('candidate_source', 'unknown')] += 1
                viral = source_candidate.get('viral_evidence')
                if viral:
                    viral_counts['candidates_with_viral_evidence'] += 1
                    if 'cytomegalovirus' in name.lower() or name.lower() == 'cmv':
                        viral_counts['cmv_candidates'] += 1
                all_blockers.update(candidate['reporting'].get('blockers') or [])
                details = {
                    'cohort': cohort, 'patient_id': pid, 'organism_name': name,
                    'tier': tier, 'candidate_source': candidate.get('candidate_source'),
                    'evidence_mode': profile, 'decision_path': candidate['decision_path'],
                    'reporting': candidate['reporting'], 'clinical_evidence': candidate['clinical_evidence'],
                    'analytical_evidence': candidate['analytical_evidence'],
                    'viral_evidence': viral, 'source_file': str(path.relative_to(root)),
                }
                if not selected and correct:
                    audit['missed_in_context'].append(details)
                if selected and not correct:
                    audit['selected_false_positives'].append(details)
            require(set(case_selected) == {c['organism_name'] for c in payload['final_report']}
                    == {c['organism_name'] for c in packet['complete_report']}, 'Final report mismatch')
            metrics_row = next(r for r in patient_labels if r['cohort'] == cohort and r['patient_id'] == pid)
            require(set(case_selected) == {s.strip() for s in metrics_row['after_integrated'].split(';') if s.strip()},
                    'Patient evaluation predictions mismatch')
            answer_names = [s.strip() for s in metrics_row['answers'].split(';') if s.strip()]
            if metrics_row['label_status'] != 'positive':
                answer_names = []
            truth = {diagnostic_name(s): s for s in answer_names}
            predicted = {diagnostic_name(s) for s in case_selected}
            pool = {diagnostic_name(c['organism_name']) for c in candidates}
            require(len(truth) == len(answer_names) and len(predicted) == len(case_selected),
                    'Diagnostic reconciliation cannot collapse distinct labels or predictions')
            fn = set(truth) - predicted
            missing_in_pool = fn & pool
            missing_outside = fn - pool
            distinct_in_pool_fn += len(missing_in_pool)
            distinct_outside_pool_fn += len(missing_outside)
            reconstructed_tp += len(set(truth) & predicted)
            reconstructed_fp += len(predicted - set(truth))
            reconstructed_fn += len(fn)
            item = {'cohort': cohort, 'patient_id': pid, 'candidate_count': len(candidates),
                    'tiers': dict(case_tiers), 'evidence_mode': profile, 'label_status': metrics_row['label_status'],
                    'provided_answers': metrics_row['answers'], 'selected': case_selected,
                    'diagnostic_fn_in_pool': [truth[k] for k in sorted(missing_in_pool)],
                    'diagnostic_fn_outside_pool': [truth[k] for k in sorted(missing_outside)]}
            cases.append(item)
            audit['patients'].append(item)
        metrics = evaluation['metrics'][cohort]['after_integrated']
        require(selected_tp == metrics['true_positive'] and total_selected == metrics['predicted'],
                f'Candidate reference flags disagree with supplied metrics: {cohort}')
        require(total_selected - selected_tp == metrics['false_positive'], 'FP count mismatch')
        diagnostic_reproduces_official_counts = (reconstructed_tp, reconstructed_fp, reconstructed_fn) == (
            metrics['true_positive'], metrics['false_positive'], metrics['false_negative'])
        require(diagnostic_reproduces_official_counts, f'Diagnostic name mapping does not reproduce baseline: {cohort}')
        no_pathogen = [r for r in cases if r['label_status'] != 'positive']
        audit['cohorts'][cohort] = {
            'patient_count': len(cases), 'candidate_count': total_candidates,
            'metrics': metrics, 'tiers': {k: dict(v) for k, v in by_tier.items()},
            'reference_positive_context_candidate_rows': context_tp,
            'diagnostic_fn_in_pool': distinct_in_pool_fn,
            'diagnostic_fn_outside_pool': distinct_outside_pool_fn,
            'diagnostic_baseline_counts_reproduce_supplied_metrics': diagnostic_reproduces_official_counts,
            'diagnostic_pool_recall_ceiling': (selected_tp + distinct_in_pool_fn) / metrics['answer_count'],
            'evidence_profiles': dict(profiles), 'no_pathogen_cases': no_pathogen,
            'candidate_free_cases': [r['patient_id'] for r in cases if not r['candidate_count']],
        }
    require(observed_keys == set(label_by_key), 'Candidate coverage mismatch')
    audit['tier_totals'] = {k: dict(v) for k, v in tier_total.items()}
    audit['candidate_source_counts'] = dict(candidate_sources)
    audit['decision_evidence_modes'] = dict(evidence_modes)
    audit['viral_counts'] = dict(viral_counts)
    audit['schema_observations'] = {'patient_json_key_presence': dict(shape_keys),
                                    'old_merge_key_presence': dict(old_merge_keys)}
    audit['reporting_blocker_counts'] = dict(all_blockers)
    audit['cross_cohort_numeric_patient_id_collisions'] = {k: v for k, v in all_case_ids.items() if len(v) > 1}
    audit['pooled_metrics'] = {
        k: sum(v['metrics'][k] for v in audit['cohorts'].values())
        for k in ['true_positive', 'false_positive', 'false_negative', 'predicted', 'answer_count']
    }
    pooled = audit['pooled_metrics']
    pooled['precision'] = pooled['true_positive'] / pooled['predicted']
    pooled['recall'] = pooled['true_positive'] / pooled['answer_count']
    pooled['f1'] = 2 * pooled['true_positive'] / (pooled['predicted'] + pooled['answer_count'])
    audit['reference_positive_context_candidate_row_count'] = len(audit['missed_in_context'])
    audit['diagnostic_fn_in_pool'] = sum(v['diagnostic_fn_in_pool'] for v in audit['cohorts'].values())
    audit['diagnostic_fn_outside_pool'] = sum(v['diagnostic_fn_outside_pool'] for v in audit['cohorts'].values())
    require(digest(ARCHIVE.read_bytes()) == archive_hash, 'Archive changed')
    output = HERE / 'release_audit.json'
    write_new(output, (json.dumps(audit, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))
    print(json.dumps({k: audit[k] for k in [
        'source_file_count', 'manifest_entries_verified', 'tier_totals', 'decision_evidence_modes',
        'pooled_metrics', 'reference_positive_context_candidate_row_count',
        'diagnostic_fn_in_pool', 'diagnostic_fn_outside_pool', 'schema_observations', 'viral_counts']},
        ensure_ascii=True, indent=2))
    print('COHORTS ' + json.dumps({k: {x: v[x] for x in ['candidate_count',
           'reference_positive_context_candidate_rows', 'diagnostic_fn_in_pool', 'diagnostic_fn_outside_pool',
           'diagnostic_pool_recall_ceiling']} for k, v in audit['cohorts'].items()}, ensure_ascii=True))


if __name__ == '__main__':
    main()
