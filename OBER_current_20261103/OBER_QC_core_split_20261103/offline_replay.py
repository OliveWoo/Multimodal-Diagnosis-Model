"""Offline-only ablation; original saved model outputs and files are read-only."""
import argparse
import copy
import csv
import io
import json
import sys
from pathlib import Path
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / 'OBER_QC_repeat_20261103'))
import repeat_driver as source
from split_policy import core_profile, split_consensus

s, runner = source.s, source.runner
BUNDLE = HERE / 'runs/offline_v1'
STRATEGY = 'QC_SPLIT_NEGATIVE'
PRIMARY = f'{STRATEGY}_F0_match_k2'
BASE_PRIMARY = 'QC_CONSENSUS_F0_match_k2'


def load(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def write_json(path, value):
    s.write_new(path, (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + '\n').encode('utf-8'))


def write_csv(path, rows):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else [])
    writer.writeheader()
    writer.writerows(rows)
    s.write_new(path, stream.getvalue().encode('utf-8-sig'))


def freeze(path, files):
    write_json(path, {'files': {str(p.resolve()): s.digest(p.read_bytes()) for p in sorted(set(files))},
                      'offline_replay_only': True, 'new_API_calls': 0, 'independent_test': False})


def source_files(manifest):
    source.verify_file_manifest(manifest)
    return [Path(p) for p in load(manifest)['files']] + [manifest]


def verify_sources():
    source.verify_file_manifest(source.BUNDLE / 'comparison_freeze.json')
    for label, folder in source.ROUNDS.items():
        for stage in ('inputs', 'QC_outputs', 'article_decisions', 'A1_inputs', 'A1_outputs', 'predictions', 'evaluation'):
            source.verify_file_manifest(folder / (stage + '_freeze.json'))


def bind(label):
    s.RUN = BUNDLE / ('round_' + label)
    s.c.RUN = s.RUN
    return s.RUN


def prepare():
    verify_sources()
    if (BUNDLE / 'setup_freeze.json').exists():
        source.verify_file_manifest(BUNDLE / 'setup_freeze.json')
        return
    for label, parent in source.ROUNDS.items():
        dest = bind(label)
        plan = copy.deepcopy(load(parent / 'plan.json'))
        plan.update(primary=PRIMARY, strategies=[STRATEGY], new_API_calls=0,
                    offline_core_split=True, experimental_not_default=True)
        records = copy.deepcopy(load(parent / 'article_decisions_without_answers.json'))
        article_files, input_files = [], []
        for article in plan['articles']:
            old = records[article['key']]
            baseline = old['strategies']['QC_CONSENSUS']
            if old['normalized_pair']:
                left, right = old['normalized_pair']
                s.require(s.consensus(left, right) == baseline, 'Original consensus not reproducible')
                derived = split_consensus(left, right, baseline)
            else:
                derived = copy.deepcopy(baseline)
                derived.update(core_recovered_from_detail_disagreement=False,
                               detail_QC_pending=False, detail_QC_fields=[])
            s.require(derived['qualified'] is not True or baseline['qualified'] is True,
                      'Split must never create a new qualified positive')
            old['strategies'][STRATEGY] = derived
            for name in ('source.json', 'routing.json'):
                new_path = dest / article['folder'] / name
                s.write_new(new_path, (parent / article['folder'] / name).read_bytes())
                input_files.append(new_path)
            path = dest / article['folder'] / 'article_decisions.json'
            write_json(path, old)
            article_files.append(path)
        write_json(dest / 'plan.json', plan)
        write_json(dest / 'A1_jobs.json', load(parent / 'A1_jobs.json'))
        write_json(dest / 'article_decisions_without_answers.json', records)
        deps = list(HERE.glob('*.py')) + [HERE / 'README.md', Path(runner.__file__), Path(s.__file__)]
        freeze(dest / 'inputs_freeze.json', input_files + deps + [dest / 'plan.json', parent / 'inputs_freeze.json'])
        freeze(dest / 'QC_outputs_freeze.json', source_files(parent / 'QC_outputs_freeze.json'))
        freeze(dest / 'article_decisions_freeze.json', article_files + [dest / 'article_decisions_without_answers.json'])
        freeze(dest / 'A1_inputs_freeze.json', source_files(parent / 'A1_inputs_freeze.json') + [dest / 'A1_jobs.json'])
        freeze(dest / 'A1_outputs_freeze.json', source_files(parent / 'A1_outputs_freeze.json'))
    freeze(BUNDLE / 'setup_freeze.json', list(HERE.glob('*.py')) + [HERE / 'README.md'] + [
        BUNDLE / ('round_' + label) / (stage + '_freeze.json')
        for label in ('A', 'B') for stage in ('inputs', 'QC_outputs', 'article_decisions', 'A1_inputs', 'A1_outputs')])


def predict():
    prepare()
    before = runner.STRATEGIES
    try:
        runner.STRATEGIES = [STRATEGY]
        # These are fail-closed guards, not mocks returning invented model output.
        with patch.object(s.c, 'infer', side_effect=RuntimeError('API forbidden in offline experiment')), \
             patch.object(s.c, 'dispatch', side_effect=RuntimeError('API dispatch forbidden')):
            for label in ('A', 'B'):
                bind(label)
                runner.predict()
    finally:
        runner.STRATEGIES = before
    freeze(BUNDLE / 'both_predictions_freeze.json', [
        BUNDLE / ('round_' + label) / 'predictions_freeze.json' for label in ('A', 'B')])
    verify_sources()
    print('Both offline predictions frozen; no new API calls and no reference labels read.', flush=True)


def selected(rows, setting):
    return {(r['case_id'], r['organism']) for r in rows if r['decisions'][setting]}


def evaluate():
    verify_sources()
    for stage in ('setup', 'both_predictions'):
        source.verify_file_manifest(BUNDLE / (stage + '_freeze.json'))
    predictions, decisions, baseline_predictions = {}, {}, {}
    for label, parent in source.ROUNDS.items():
        dest = bind(label)
        for stage in ('inputs', 'QC_outputs', 'article_decisions', 'A1_inputs', 'A1_outputs', 'predictions'):
            source.verify_file_manifest(dest / (stage + '_freeze.json'))
        predictions[label] = load(dest / 'predictions_without_answers.json')
        decisions[label] = load(dest / 'article_decisions_without_answers.json')
        baseline_predictions[label] = load(parent / 'predictions_without_answers.json')

    # First reference read in this ablation, after both rounds' predictions freeze.
    from inspect_release import diagnostic_name as canon, POSITIVE_TIERS
    reference = ROOT / 'OBER_v7_review_20261030/cohort49_review.json'
    refs = load(reference)['cases']
    truth = {r['case_id']: {canon(n) for n in r['reference']} for r in refs}
    groups = {r['case_id']: r['analysis_group'] for r in refs}
    clinical = load(s.c.old.d.SOURCE / 'candidates_not_model_input.json')
    base = {cid: {canon(r['organism']) for r in clinical if r['case_id'] == cid and r['baseline_tier'] in POSITIVE_TIERS}
            for cid in truth}
    s.require(len(truth) == 49 and sum(map(len, truth.values())) == 87 and sum(map(len, base.values())) == 97,
              'Changed cohort, reference or original formal predictions')
    metric = s.c.old.d.clinical.metrics

    def score(pairs, group):
        tp = fp = fn = 0
        for cid, answer in truth.items():
            if group != 'ALL49' and groups[cid] != group:
                continue
            predicted = base[cid] | {canon(n) for c, n in pairs if c == cid}
            tp += len(predicted & answer)
            fp += len(predicted - answer)
            fn += len(answer - predicted)
        return metric(tp, fp, fn)

    metrics, curves, candidates, article_rows, coverage = [], [], [], [], {}
    for label in ('A', 'B'):
        rows = predictions[label]
        base_rows = baseline_predictions[label]
        baseline_score = score(selected(base_rows, BASE_PRIMARY), 'ALL49')
        expected = load(source.ROUNDS[label] / 'evaluation.json')['primary']
        for field in ('TP', 'FP', 'FN', 'precision', 'recall', 'F1'):
            s.require(abs(baseline_score[field] - expected[field]) < 1e-12, 'Baseline scorer mismatch')
        for setting in rows[0]['decisions']:
            for group in ('ALL49', 'KH33', 'KH5', 'A11'):
                metrics.append(dict(round=label, setting=setting, group=group,
                                    **score(selected(rows, setting), group)))
        for k in range(11):
            old_setting, new_setting = f'QC_CONSENSUS_F0_match_k{k}', f'{STRATEGY}_F0_match_k{k}'
            old_picks, new_picks = selected(base_rows, old_setting), selected(rows, new_setting)
            curves.append(dict(round=label, k=k, original_F1=score(old_picks, 'ALL49')['F1'],
                               split_F1=score(new_picks, 'ALL49')['F1'],
                               changed_selection_pairs=len(old_picks ^ new_picks)))
        oldmap = {r['token']: r for r in base_rows}
        for row in rows:
            if not row.get('variants'):
                continue
            oldrow = oldmap[row['token']]
            ov, nv = oldrow['variants']['QC_CONSENSUS_F0'], row['variants'][f'{STRATEGY}_F0']
            out = dict(round=label, case_id=row['case_id'], organism=row['organism'],
                       original_selected=oldrow['decisions'][BASE_PRIMARY], split_selected=row['decisions'][PRIMARY],
                       original_S=ov['S_pass'], split_S=nv['S_pass'])
            for arm in ('G', 'T'):
                def count(value):
                    c = value['counts']
                    return None if c is None else c[arm]['strong_match_count'] + c[arm]['partial_match_count']
                out[f'original_{arm}_match'], out[f'split_{arm}_match'] = count(ov), count(nv)
            candidates.append(out)
        items = list(decisions[label].values())
        olditems = [r['strategies']['QC_CONSENSUS'] for r in items]
        newitems = [r['strategies'][STRATEGY] for r in items]
        coverage[label] = dict(
            original_core_pending=sum(r['status'] == 'unresolved_qc' for r in olditems),
            split_core_pending=sum(r['status'] == 'unresolved_qc' for r in newitems),
            recovered_negative_core=sum(r['core_recovered_from_detail_disagreement'] for r in newitems),
            detail_QC_pending=sum(r['detail_QC_pending'] for r in newitems),
            original_qualified=sum(r['qualified'] is True for r in olditems),
            split_qualified=sum(r['qualified'] is True for r in newitems))
        for entry in items:
            old, new = entry['strategies']['QC_CONSENSUS'], entry['strategies'][STRATEGY]
            article_rows.append(dict(round=label, pmid=entry['pmid'], organism=entry['organism'],
                                     original_status=old['status'], split_status=new['status'],
                                     original_qualified=old['qualified'], split_qualified=new['qualified'],
                                     recovered_negative_core=new['core_recovered_from_detail_disagreement'],
                                     detail_QC_pending=new['detail_QC_pending'],
                                     detail_fields='|'.join(new['detail_QC_fields'])))

    old_flips = new_flips = 0
    old_differences = new_differences = 0
    for key, entry in decisions['A'].items():
        other = decisions['B'][key]
        for strategy in ('QC_CONSENSUS', STRATEGY):
            a, b = entry['strategies'][strategy], other['strategies'][strategy]
            flip = (a['status'] == 'unresolved_qc') != (b['status'] == 'unresolved_qc')
            different = core_profile(a) != core_profile(b)
            if strategy == STRATEGY:
                new_flips += flip
                new_differences += different
            else:
                old_flips += flip
                old_differences += different
    primary = {label: score(selected(predictions[label], PRIMARY), 'ALL49') for label in ('A', 'B')}
    diffs = {label: sorted(selected(predictions[label], PRIMARY) ^ selected(baseline_predictions[label], BASE_PRIMARY))
             for label in ('A', 'B')}
    result = dict(primary=primary, coverage=coverage, changed_primary_pairs_from_baseline=diffs,
                  same_primary_selection_between_rounds=selected(predictions['A'], PRIMARY) == selected(predictions['B'], PRIMARY),
                  cross_round_core_status_flips_before=old_flips, cross_round_core_status_flips_after=new_flips,
                  cross_round_core_profile_differences_before=old_differences,
                  cross_round_core_profile_differences_after=new_differences,
                  selected={label: sorted(selected(predictions[label], PRIMARY)) for label in ('A', 'B')},
                  new_API_calls=0, source_model_outputs_unchanged=True, original_formal_pairs_preserved=97,
                  reference_sha256=s.digest(reference.read_bytes()), independent_test=False,
                  experimental_not_adopted_as_default=True, primary_setting=PRIMARY)
    write_json(BUNDLE / 'comparison.json', result)
    for name, values in [('metrics', metrics), ('k_curve', curves), ('candidate_changes', candidates), ('article_comparison', article_rows)]:
        write_csv(BUNDLE / (name + '.csv'), values)
    report = ['# 核心／細節QC分離：保守離線試驗結果', '',
              '原版本維持使用者目前可接受的基準；本實驗不自動採用。使用既有A/B兩輪輸出，新增API呼叫0次。', '',
              '## 固定主設定k2', '', '|輪次|版本|TP|FP|FN|Precision|Recall|F1|核心未解QC|細節QC|',
              '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for label in ('A', 'B'):
        for version, m, pending, detail in [
            ('原版', score(selected(baseline_predictions[label], BASE_PRIMARY), 'ALL49'), coverage[label]['original_core_pending'], '未獨立拆分'),
            ('保守拆分', primary[label], coverage[label]['split_core_pending'], coverage[label]['detail_QC_pending'])]:
            report.append(f"|{label}|{version}|{m['TP']}|{m['FP']}|{m['FN']}|{m['precision']:.2%}|{m['recall']:.2%}|{m['F1']:.2%}|{pending}|{detail}|")
    report += ['', f'- A輪恢復{coverage["A"]["recovered_negative_core"]}篇負向核心結論；B輪恢復{coverage["B"]["recovered_negative_core"]}篇。細節分歧仍保留，並非被消除。',
               f'- 跨輪「核心已解／未解」翻轉：{old_flips} → {new_flips}篇。',
               f'- 跨輪核心五欄／未解狀態差異：{old_differences} → {new_differences}篇。這不是拿原全欄位分類差異7篇直接比較。',
               f'- 每輪相對原k2清單的差異：A {len(diffs["A"])}筆；B {len(diffs["B"])}筆。兩輪新清單相同：{result["same_primary_selection_between_rounds"]}。', '',
               '## k敏感度（不換主k2）', '', '|k|原A F1|拆分A F1|原B F1|拆分B F1|', '|---:|---:|---:|---:|---:|']
    for k in range(11):
        a = next(r for r in curves if r['round'] == 'A' and r['k'] == k)
        b = next(r for r in curves if r['round'] == 'B' and r['k'] == k)
        report.append(f"|{k}|{a['original_F1']:.2%}|{a['split_F1']:.2%}|{b['original_F1']:.2%}|{b['split_F1']:.2%}|")
    report += ['', '## 限制與保留決策', '',
               '- 未重新呼叫LLM：這是既有答案的確定性彙整規則重播，不是兩次新的穩定性驗證。',
               '- 原58篇來源、G/T、A1病例六欄、固定答案、原97筆正式選擇均不變；核心恢復可能改變A1衍生文獻欄位，因此有完整重播而非假設下游不變。',
               '- 只有雙方明確一致不合格的核心可恢復；正向支持與核心分歧繼續待查，沒有新的正向文獻資格。',
               '- 細節QC另列，不可把核心待查減少當作模型分類正確或不確定性消失。沒有修正原模型分類本身，也不證明醫學判讀正確。',
               '- 已反覆開發的49例、87答案，不是獨立test；不更換主k、不加入病原特例。',
               '- 使用者目前可接受原結果，保留原版，不自動把小實驗當正式改版。', '',
               '使用者決策紀錄：專案根目錄OBER_CURRENT_STATUS.md。逐篇／逐候選資料見runs/offline_v1。']
    report_path = HERE / 'RESULTS_zh.md'
    s.write_new(report_path, ('\n'.join(report) + '\n').encode('utf-8'))
    freeze(BUNDLE / 'evaluation_freeze.json', [report_path, BUNDLE / 'comparison.json'] + list(BUNDLE.glob('*.csv')))
    verify_sources()
    print(json.dumps(result, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=('predict', 'evaluate'))
    {'predict': predict, 'evaluate': evaluate}[parser.parse_args().stage]()
