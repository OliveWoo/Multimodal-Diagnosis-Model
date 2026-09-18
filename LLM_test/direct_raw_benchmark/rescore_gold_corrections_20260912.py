"""Apply explicit PDF Gold corrections and rescore frozen predictions, offline.

Preserves original Gold, predictions, scoring code, prior reports, and all API
settings. Writes a versioned KH Gold and a new report directory only.
"""
from __future__ import annotations

import copy
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

import direct_raw_runner as r
import evaluate_direct_raw as ev
from input_identity import sha256_file

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
MANIFEST = BASE / "gold_corrections_20260912.json"
PRIOR = BASE / "reports/astra_20260906/comparison.json"
OLD_KH = ROOT / "RAG_re/gold/KMUH_clinical_pathogen_gold_20260811.json"
NEW_KH = ROOT / "RAG_re/gold/KMUH_clinical_pathogen_gold_corrected_20260912.json"
REPORT = BASE / "reports/gold_corrected_20260912"
COHORTS = [("all_41", "高醫＋two-hospital A：41 位"),
           ("kmuh_answer_positive_30", "高醫：30 位"),
           ("two_hospital_patient_quality_A_11", "Two-hospital A：11 位")]
MODELS = ("sol", "luna", "terra", "astra")


def same_set(left, right):
    return set(left) == set(right) and len(left) == len(set(left)) and len(right) == len(set(right))


def apply_corrections(original, manifest):
    revised = copy.deepcopy(original)
    cases = {str(c['patient_id']): c for c in revised['cases']}
    if len(cases) != len(revised['cases']):
        raise ValueError("Duplicate Gold patient")
    seen = set()
    journal = []
    for correction in manifest['corrections']:
        pid = correction['patient_id']
        if pid in seen or pid not in cases:
            raise ValueError(f"Duplicate/unknown correction patient P{pid}")
        seen.add(pid)
        case = cases[pid]
        if not same_set(case['pathogens'], correction['expected_before']):
            raise ValueError(f"Original Gold differs from reviewed version: P{pid}")
        if correction['action'] == 'confirm_unchanged':
            if not same_set(case['pathogens'], correction['after']):
                raise ValueError("Confirmation cannot change labels")
        elif correction['action'] != 'replace':
            raise ValueError("Unsupported correction action")
        if len(correction['after']) != len(set(correction['after'])):
            raise ValueError("Duplicate corrected label")
        case['pathogens'] = list(correction['after'])
        case['answer_revision'] = {
            'pdf_page': correction['pdf_page'], 'review_document_date': manifest['review_document_date'],
            'action': correction['action'], 'reason_zh': correction['reason_zh'],
            'previous_pathogens': list(correction['expected_before']),
            'raw_diagnosis_cells_semantics': 'Original workbook cells retained for provenance; not the corrected answers.',
        }
        journal.append(copy.deepcopy(correction))
    revised['development_only'] = True
    revised['answer_revision'] = {
        'id': 'PDF_20260912_applied_20260914', 'source_pdf': str((BASE / manifest['source_pdf']).resolve()),
        'source_pdf_sha256': manifest['source_pdf_sha256'], 'source_manifest': str(MANIFEST),
        'original_gold_file': str(OLD_KH), 'original_gold_sha256': sha256_file(OLD_KH),
        'applied_date': manifest['applied_date'], 'revised_patient_ids': ['8', '14', '16', '29'],
        'confirmed_unchanged_patient_ids': ['30', '39'],
        'evaluation_status': 'Post-prediction reviewed labels; not an untouched or blinded test set.',
        'historical_split_note': 'Original case split fields are retained as historical metadata, not a claim of untouched holdout validity.',
    }
    return revised, journal


def score_row(pid, cohort, gold, predicted):
    tp, fp, fn = gold & predicted, predicted - gold, gold - predicted
    precision = len(tp) / len(predicted) if predicted else 0.0
    recall = len(tp) / len(gold) if gold else 0.0
    return {'patient_id': pid, 'cohort': cohort, 'gold': sorted(gold),
            'predicted': sorted(predicted), 'tp': sorted(tp), 'fp': sorted(fp), 'fn': sorted(fn),
            'precision': precision, 'recall': recall, 'f1': ev.f1(precision, recall)}


def summarize(rows):
    return {key: ev.summarize([x for x in rows if key == 'all_41' or x['cohort'] == key])
            for key, _ in COHORTS}


def write_frozen_json(path, value):
    if path.exists():
        if r.load_json(path) != value:
            raise ValueError(f'Refusing to overwrite a different frozen artifact: {path}')
    else:
        r.write_json(path, value)


def main():
    revision = r.load_json(MANIFEST)
    pdf = (BASE / revision['source_pdf']).resolve()
    if sha256_file(pdf) != revision['source_pdf_sha256']:
        raise ValueError('Correction PDF changed')
    prior = r.load_json(PRIOR)
    old_astra = r.load_json(BASE / 'reports/astra_20260906/astra.json')
    gold_sources = [Path(p) for p in old_astra['gold_sources']]
    for path in gold_sources:
        if sha256_file(path) != prior['source_hashes']['gold'][path.name]:
            raise ValueError(f'Original Gold changed: {path.name}')
    if sha256_file(BASE / 'evaluate_direct_raw.py') != prior['source_hashes']['evaluator']:
        raise ValueError('Frozen evaluator changed')
    _, specs = r.load_manifest(r.DEFAULT_MANIFEST)
    if len(specs) != 41:
        raise ValueError('Expected exactly 41 frozen patients')
    spec_by_id = {s.numeric_id: s for s in specs}
    for correction in revision['corrections']:
        spec = spec_by_id[correction['patient_id']]
        if spec.cohort != revision['cohort'] or spec.expected_case_code != correction['case_code']:
            raise ValueError('Correction ID namespace mismatch')
    input_audits = {s.patient_id: r.build_patient_payload(s)[1] for s in specs}
    prompt_hash = r.sha256_text(r.DEFAULT_PROMPT.read_text(encoding='utf-8-sig').strip())
    if prompt_hash != prior['source_hashes']['prompt']:
        raise ValueError('Prompt changed')
    verifications = {(x['model'], x['patient_id']): x for x in prior['verified_predictions']}
    frozen_files = {}
    for name in MODELS:
        model = 'gpt-6-astra' if name == 'astra' else f'gpt-5.6-{name}'
        run = r.DEFAULT_OUTPUT_ROOT / f'direct_raw_41_{name}_idfixed_20260906'
        if {p.stem.removesuffix('_direct_raw') for p in run.glob('P*_direct_raw.json')} != {s.patient_id for s in specs}:
            raise ValueError(f'{name}: missing/extra patient output')
        for spec in specs:
            path = run / f'{spec.patient_id}_direct_raw.json'
            record = r.load_json(path)
            r.validate_cached_output(record, spec, input_audits[spec.patient_id], prompt_hash, model, 'medium', 4000)
            expected = verifications[(model, spec.patient_id)]['output_sha256']
            if sha256_file(path) != expected:
                raise ValueError(f'Frozen prediction changed: {name} {spec.patient_id}')
            frozen_files[str(path)] = expected
    old_ober_audit = r.load_json(BASE / 'reports/idfixed_20260906/comparison.json')
    ober_hashes = {x['patient_id']: x['source_sha256'] for x in old_ober_audit['OBER_per_patient']}
    decisions = {}
    for spec in specs:
        cohort = 'KH' if spec.cohort.startswith('kmuh') else 'two_hospital'
        path = ROOT / f'OBER_patient_delivery/outputs/20260905_r5_final_019feb4d/patients/{cohort}/{spec.patient_id}/evidence/final_decision.json'
        if sha256_file(path) != ober_hashes[spec.numeric_id]:
            raise ValueError(f'OBER decision changed: {spec.patient_id}')
        decisions[spec.numeric_id] = r.load_json(path)
        frozen_files[str(path)] = sha256_file(path)
    revised, journal = apply_corrections(r.load_json(OLD_KH), revision)
    documents = [revised] + [r.load_json(p) for p in gold_sources if p != OLD_KH]
    aliases = ev.build_aliases(documents)
    if aliases != ev.build_aliases([r.load_json(p) for p in gold_sources]):
        raise ValueError('Species alias rules changed')
    changed = {x['patient_id'] for x in journal if x['action'] == 'replace'}
    if changed != {'8', '14', '16', '29'}:
        raise ValueError('Unexpected correction scope')
    # All source/prediction checks finish before any versioned artifact is written.
    write_frozen_json(NEW_KH, revised)
    REPORT.mkdir(parents=True, exist_ok=True)
    write_frozen_json(REPORT / 'correction_manifest.json', revision)
    results = {}
    before_results = {}
    for name in MODELS:
        run = r.DEFAULT_OUTPUT_ROOT / f'direct_raw_41_{name}_idfixed_20260906'
        cmd = [sys.executable, str(BASE / 'evaluate_direct_raw.py'), '--run-dir', str(run),
               '--kmuh-gold', str(NEW_KH), '--out-json', str(REPORT / f'{name}.json'),
               '--out-md', str(REPORT / f'{name}.md')]
        completed = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8')
        if completed.returncode:
            raise RuntimeError(completed.stderr or completed.stdout)
        after = r.load_json(REPORT / f'{name}.json')
        before = r.load_json(BASE / f'reports/astra_20260906/{name}.json')
        previous = {x['patient_id']: x for x in before['per_patient']}
        for row in after['per_patient']:
            old = previous[row['patient_id']]
            if row['predicted'] != old['predicted']:
                raise ValueError('Prediction or matching changed')
            if row['patient_id'] not in changed and row != old:
                raise ValueError('Unreviewed patient changed')
            independently = score_row(row['patient_id'], row['cohort'], set(row['gold']), set(row['predicted']))
            if any(independently[k] != row[k] for k in independently):
                raise ValueError('Independent scoring cross-check failed')
        results[name] = after
        before_results[name] = before
    gold_by_id = {x['patient_id']: set(x['gold']) for x in results['astra']['per_patient']}
    old_gold_by_id = {x['patient_id']: set(x['gold']) for x in before_results['astra']['per_patient']}
    ober_rows, old_ober_rows = [], []
    for spec in specs:
        pid = spec.numeric_id
        predicted = {ev.canonical(x['organism'], aliases) for x in decisions[pid]['selected']}
        ober_rows.append(score_row(pid, spec.cohort, gold_by_id[pid], predicted))
        old_ober_rows.append(score_row(pid, spec.cohort, old_gold_by_id[pid], predicted))
    old_ober_metrics = summarize(old_ober_rows)
    if old_ober_metrics != prior['metrics']['上游＋OBER R5（原決策）']:
        raise ValueError('Original OBER metrics not reproduced')
    results['ober'] = {'schema_version': 'gold_correction.ober_evaluation.v1',
                       'metrics': summarize(ober_rows), 'per_patient': ober_rows,
                       'gold_sources': results['astra']['gold_sources'],
                       'prediction_source': 'Unmodified frozen upstream + OBER R5 final decisions'}
    before_results['ober'] = {'metrics': old_ober_metrics, 'per_patient': old_ober_rows}
    r.write_json(REPORT / 'ober.json', results['ober'])
    all_methods = ('ober',) + MODELS
    comparison = {}
    for name in all_methods:
        comparison[name] = {'before': before_results[name]['metrics'], 'after': results[name]['metrics']}
        if results[name]['metrics']['two_hospital_patient_quality_A_11'] != before_results[name]['metrics']['two_hospital_patient_quality_A_11']:
            raise ValueError('Two-hospital A changed despite no corrections')
        if results[name]['metrics']['all_41']['gold_pair_count'] != 76:
            raise ValueError('Expected unchanged total of 76 Gold pairs')
    changes = []
    for name in all_methods:
        previous = {x['patient_id']: x for x in before_results[name]['per_patient']}
        for row in results[name]['per_patient']:
            if row['patient_id'] in changed:
                old = previous[row['patient_id']]
                changes.append({'method': name, 'patient_id': row['patient_id'], 'before_gold': old['gold'],
                                'after_gold': row['gold'], 'predicted': row['predicted'],
                                'before': {k: old[k] for k in ('tp', 'fp', 'fn')},
                                'after': {k: row[k] for k in ('tp', 'fp', 'fn')}})
    for path, expected in frozen_files.items():
        if sha256_file(Path(path)) != expected:
            raise ValueError('Prediction changed during scoring')
    for path in gold_sources:
        if sha256_file(path) != prior['source_hashes']['gold'][path.name]:
            raise ValueError('Original Gold modified during scoring')
    report = {'schema_version': 'direct_raw_benchmark.pdf_gold_revision_comparison.v1',
              'created_at_utc': datetime.now(timezone.utc).isoformat(), 'source_pdf': str(pdf),
              'source_pdf_sha256': sha256_file(pdf), 'corrected_kmuh_gold': str(NEW_KH),
              'corrected_kmuh_gold_sha256': sha256_file(NEW_KH), 'corrections': journal,
              'methods': comparison, 'per_patient_changes': changes,
              'unchanged_prediction_file_count': len(frozen_files),
              'unchanged_prediction_sha256': frozen_files,
              'matching': 'unchanged exact species plus frozen synonyms',
              'notes': ['Gold-only post-prediction re-adjudication; no API requests or model reruns.',
                        'The PDF includes model outputs; this is not blinded independent validation.',
                        'P30 remains included as instructed despite a non-pulmonary rationale in the PDF.',
                        'P8 and P39 reads discrepancies are not changes to Gold and were not applied to raw input.',
                        'Same 30 KH + 11 two-hospital A patients; broad Aspergillosis handling unchanged.']}
    r.write_json(REPORT / 'comparison.json', report)
    titles = {'ober': '上游＋OBER R5', 'sol': 'GPT Sol', 'luna': 'GPT Luna', 'terra': 'GPT Terra', 'astra': 'GPT Astra'}
    lines = ['# 2026-09-12 PDF 更正答案後重評分', '',
             '依據 2026-9-12.pdf 六頁的明確更正／確認，僅更新高醫 Gold；不重跑模型、不改動選菌、候選池、門檻、raw reads 或同義詞。', '',
             '## 答案更正', '', '| 病例 | PDF 頁 | 原答案 | 更正後答案 | 處理 |', '|---|---:|---|---|---|']
    for item in journal:
        lines.append(f"| P{item['patient_id']} | {item['pdf_page']} | {'；'.join(item['expected_before'])} | {'；'.join(item['after'])} | {'更正' if item['action']=='replace' else '確認保留、不變'} |")
    lines += ['', '四位答案有變、兩位確認不變。新增 3 個標籤、移除 3 個標籤，因此高醫仍為 58 個 Gold 配對，two-hospital A 仍為 18 個，合計 76 個。', '']
    for cohort, title in COHORTS:
        lines += [f'## {title}', '', '| 方法 | TP | FP | FN | Recall | F1 | Precision |', '|---|---:|---:|---:|---:|---:|---:|']
        for name in all_methods:
            m = comparison[name]['after'][cohort]
            lines.append(f"| {titles[name]} | {m['tp']} | {m['fp']} | {m['fn']} | {ev.percent(m['recall'])} | {ev.percent(m['f1'])} | {ev.percent(m['precision'])} |")
        lines += ['']
    lines += ['## 高醫＋A：更正前 → 更正後', '', '| 方法 | Recall | F1 | Precision |', '|---|---:|---:|---:|']
    for name in all_methods:
        b, a = comparison[name]['before']['all_41'], comparison[name]['after']['all_41']
        lines.append(f"| {titles[name]} | {ev.percent(b['recall'])} → {ev.percent(a['recall'])} | {ev.percent(b['f1'])} → {ev.percent(a['f1'])} | {ev.percent(b['precision'])} → {ev.percent(a['precision'])} |")
    lines += ['', '## 核對與限制', '',
              '- 評分單位是 (病人, 病原) 配對的 pooled micro，不把不同病人的同名菌合併為一個。',
              '- Candida albicans、Candida tropicalis、Candida dubliniensis 仍為不同菌種；其餘 frozen synonyms 與 Aspergillosis 等廣義標籤處理皆不變。',
              '- 164 份 GPT 回答＋41 份 OBER 最終決策雜湊驗證不變；三份原始 Gold 檔案、舊報告與舊模型輸出均保留。',
              '- PDF P39 的 E. faecium reads=83 只是資料差異說明，沒有新增該菌為正確答案；P8 的 reads 差異也未回填到任何輸入。',
              '- PDF P30 的保留理由涉及腹內感染／造口傷口；此處依文件保留 B. fragilis，不代表已證明它是肺部來源。',
              '- PDF 含原模型結果，因此更正屬於看過預測後的答案再審查，不是盲審或 untouched test。分數變化來自 Gold 更正，不是模型重新訓練或進步。',
              '- 不處理本次 PDF 以外的資料疑慮、病人分級或標籤爭議。', '',
              '## 檔案與重現', '',
              '- [更正後高醫 Gold](../../../../RAG_re/gold/KMUH_clinical_pathogen_gold_corrected_20260912.json)',
              '- [完整比較與逐病人變化](comparison.json)',
              '- [OBER 逐病人結果](ober.json)',
              '- [Sol](sol.md)／[Luna](luna.md)／[Terra](terra.md)／[Astra](astra.md)', '',
              '在 LLM_test/direct_raw_benchmark 執行 `python rescore_gold_corrections_20260912.py` 可離線重算，不會呼叫 API。', '']
    (REPORT / 'comparison_zh.md').write_text('\n'.join(lines), encoding='utf-8')
    print(json.dumps({'corrected_gold': str(NEW_KH), 'metrics': {titles[n]: comparison[n]['after']['all_41'] for n in all_methods}, 'report': str(REPORT / 'comparison_zh.md')}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
