"""A1 article extraction -> unchanged legacy truth table -> deterministic rescue.

No candidate-level LLM recommendation is requested or accepted. prepare/run do
not read reference answers. All decisions across k are frozen before evaluate.
"""
import argparse
import concurrent.futures
import copy
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

import information_pilot_recovery as information
from api_probe import environment
from inspect_release import HERE, digest, load, require, write_new

ROOT = HERE.parent
RUN = HERE / 'runs/pilot_v6_A1_rules_20261031'
LEGACY = ROOT / 'RAG_re_casefit'
sys.path.insert(0, str(LEGACY))
from rag_re_casefit.rules import derive_article_verdict, apply_traceability_gate, aggregate_judgments

old = information.p.old
v4 = information.p.v4
MODEL = old.MODEL
FIELDS = ['organism_scope', 'human_clinical_context', 'site_alignment', 'syndrome_alignment',
          'specimen_alignment', 'host_alignment', 'phenotype_alignment', 'temporal_alignment']
ALIGNMENTS = ['match', 'partial', 'mismatch', 'unknown']
K_VALUES = list(range(1, 11))
B_RANK_LIMIT = 3
ACCEPTED_QC = ['evaluable', 'partial_evaluable']
SOURCE_PROMPT = LEGACY / 'prompts/article_casefit_v1.txt'
PROMPT = SOURCE_PROMPT.read_text(encoding='utf-8').strip() + '''

The original categorical definitions above are unchanged. This input provides
source patient facts rather than inferred host tags. The protocol target is a
study question, not proof the patient has pneumonia. Disease site and syndrome
alignment can cite P_SITE/P_SYNDROME as the comparison target, explicitly noting
their protocol role. All other alignment claims need actual patient fact IDs.
The index specimen is a sampling fact, not proof of the candidate's causality.
An image marked timing_unknown cannot by itself establish this episode's phenotype
or timing: use unknown if no episode-linked source establishes that feature.
Recorded steroid exposure does not establish prolonged treatment unless stated.
Do not infer neutropenia from a malignancy label. Preserve contradictory records.
Patient microbiology results, mNGS strength, upstream decisions, and answers are
withheld. Do not infer them. This is not a request for diagnosis or selection.
For every categorical field return separate article sentence IDs and patient fact
IDs. organism_scope/human_clinical_context need article evidence only; a known
alignment (match/partial/mismatch) needs BOTH sides. Unknown may have empty IDs.
Use only this article and these facts. Exact article excerpts remain <=20 words.
Return the requested extraction schema only, never an admission recommendation,
final verdict, probability of causality, or a private chain of thought.
Keep the short rationale in Traditional Chinese.
'''


def obj(properties):
    return old.obj(properties)


def ids():
    return {'type': 'array', 'items': {'type': 'string'}}


SCHEMA = obj({
    'candidate_token': {'type': 'string'}, 'pmid': {'type': 'string'},
    'organism_scope': {'type': 'string', 'enum': ['exact_species', 'genus_only', 'other', 'unclear']},
    'human_clinical_context': {'type': 'string', 'enum': ['causal_infection', 'colonization_or_contamination',
        'reactivation_or_bystander', 'detection_only', 'nonhuman', 'unclear']},
    **{key: {'type': 'string', 'enum': ALIGNMENTS} for key in FIELDS[2:]},
    'field_evidence': obj({key: obj({'article_sentence_ids': ids(), 'patient_fact_ids': ids()}) for key in FIELDS}),
    'support_span': {'type': 'string'}, 'counter_span': {'type': 'string'},
    'rationale': {'type': 'string'},
})
CLINICAL_NAMES = {
    'Human alphaherpesvirus 1': ['HSV-1', 'human herpesvirus 1', 'herpes simplex virus type 1'],
    'Human betaherpesvirus 5': ['human cytomegalovirus', 'CMV'],
}


def save(path, value):
    write_new(path, old.encoded(value))


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def b_strength(branches):
    """Relative analytical strength only. No fake mapping to historical B_STRICT."""
    checks = []
    for i, row in enumerate(branches):
        rank = row.get('rank_in_retained_universe'); reads = row.get('reads')
        qc = row.get('qc_status')
        if finite(reads) and reads <= 0:
            state = False
        elif qc and qc not in ACCEPTED_QC:
            state = False
        elif finite(rank) and rank > B_RANK_LIMIT:
            state = False
        elif finite(rank) and rank >= 1 and finite(reads) and reads > 0 and qc in ACCEPTED_QC:
            state = True
        else:
            state = None
        checks.append({'branch_index': i, 'seq_id': row.get('seq_id'), 'reads': reads,
            'rank_in_retained_universe': rank, 'analytical_rank_priority_not_used_as_rank': row.get('analytical_rank_priority'),
            'rpm_total_not_thresholded': row.get('rpm_total'), 'qc_status': qc, 'passes': state})
    positive = True if any(c['passes'] is True for c in checks) else None if not checks or any(c['passes'] is None for c in checks) else False
    return {'version': 'B_v7_relative_rank3', 'positive': positive, 'rank_limit': B_RANK_LIMIT,
        'branches': checks, 'old_B_STRICT_reproducible': False,
        'old_missing_fields': ['mngs_signal_tier', 'reads_tier', 'dominance_tier', 'reads_percentile'],
        'interpretation': 'At least one SAME branch has positive reads, accepted QC and retained-universe rank <=3. Relative analytical research gate only; no cross-branch sum, independent confirmation or calibrated diagnostic cutoff.'}


def a1_view(core):
    # Reuse verified sources; no new LLM case summary, no test-result leakage into A.
    facts = [copy.deepcopy(f) for f in core['facts'] if f['kind'] in ['index', 'background', 'imaging']]
    facts += [{'id': 'P_SITE', 'kind': 'protocol_target', 'temporal_role': 'protocol',
               'value': 'lower_respiratory', 'caveat': 'Study target, not patient-level confirmation.'},
              {'id': 'P_SYNDROME', 'kind': 'protocol_target', 'temporal_role': 'protocol',
               'value': 'lower_respiratory_infection', 'caveat': 'Study target, not patient-level confirmation.'}]
    return {'candidate_token': core['candidate_token'], 'target_syndrome': 'lower_respiratory_infection',
        'clinical_site': 'lower_respiratory', 'index_specimen_type': core['index_specimen_type'],
        'time_semantics': core['time_semantics'], 'absence_semantics': core['absence_semantics'], 'facts': facts,
        'withheld': ['mNGS reads/RPM/rank', 'culture/PCR/FilmArray results', 'upstream selection and rationale', 'reference answers']}


def article_with_sentences(article):
    return old.with_sentences([article])[0]


def validate_extraction(output, view, article):
    require(set(output) == set(SCHEMA['properties']), 'Unexpected output fields; candidate recommendations are not allowed')
    require(set(output['field_evidence']) == set(FIELDS), 'Missing per-field evidence mapping')
    require(output['candidate_token'] == view['candidate_token'] and output['pmid'] == article['pmid'], 'Wrong article/candidate binding')
    facts = {f['id']: f for f in view['facts']}
    sentences = {s['id']: s['text'] for s in article['sentences']}
    failures = []
    for key in FIELDS:
        refs = output['field_evidence'][key]
        require(set(refs) == {'article_sentence_ids', 'patient_fact_ids'}, 'Unexpected evidence fields')
        a, f = refs['article_sentence_ids'], refs['patient_fact_ids']
        if not set(a) <= set(sentences) or not set(f) <= set(facts):
            failures.append(key + ':unknown_source_id'); continue
        known = output[key] not in ['unknown', 'unclear']
        if known and not a:
            failures.append(key + ':missing_article_evidence')
        if key.endswith('_alignment') and known:
            if not f:
                failures.append(key + ':missing_patient_evidence')
            if key not in ['site_alignment', 'syndrome_alignment'] and f and all(facts[i]['kind'] == 'protocol_target' for i in f):
                failures.append(key + ':protocol_is_not_patient_fact')
            if key in ['phenotype_alignment', 'temporal_alignment'] and f and all(
                    facts[i].get('temporal_role') == 'timing_unknown' for i in f):
                failures.append(key + ':only_unresolved_timing_context')
    blob = ' '.join((article['title'] + ' ' + article['abstract']).split())
    def span_valid(text):
        compact = ' '.join(text.split())
        return bool(compact and len(compact.split()) <= 20 and compact.casefold() in blob.casefold())
    extraction = {key: output[key] for key in FIELDS}
    decision = derive_article_verdict(extraction)
    support_valid, counter_valid = span_valid(output['support_span']), span_valid(output['counter_span'])
    decision = apply_traceability_gate(decision, support_span_valid=support_valid, counter_span_valid=counter_valid)
    if failures:
        decision.update({'verdict': 'insufficient', 'counterevidence': False, 'retrieval_mismatch': False,
                         'traceability_status': 'invalid_field_citations'})
    return {'extraction': output, **decision, 'support_span_valid': support_valid,
            'counter_span_valid': counter_valid, 'field_trace_failures': failures,
            'traceability_valid': not failures and decision['traceability_status'] == 'valid'}


def no_abstract(view, article):
    output = {'candidate_token': view['candidate_token'], 'pmid': article['pmid'],
        'organism_scope': 'unclear', 'human_clinical_context': 'unclear',
        **{key: 'unknown' for key in FIELDS[2:]},
        'field_evidence': {key: {'article_sentence_ids': [], 'patient_fact_ids': []} for key in FIELDS},
        'support_span': '', 'counter_span': '', 'rationale': '沒有可用摘要；保留不足，不以標題補造病例條件。'}
    return {'status': 'no_abstract', 'pmid': article['pmid'], **validate_extraction(output, view, article), 'usage': None}


def rule_configs():
    rules = [{'id': 'E0', 'family': 'baseline', 'count': None, 'k': None},
             {'id': 'B_only', 'family': 'B', 'count': None, 'k': None}]
    for count in ['strong', 'match']:
        for k in K_VALUES:
            for family in ['A', 'AB']:
                rules.append({'id': f'{family}_{count}_k{k}', 'family': family, 'count': count, 'k': k})
    return rules


def admit(row, rule):
    if row['baseline_tier'] != 'Context' or rule['family'] == 'baseline':
        return False
    if rule['family'] == 'B':
        return row['B']['positive'] is True
    count = row['A1']['strong_match_count']
    if rule['count'] == 'match':
        count += row['A1']['partial_match_count']
    return count >= rule['k'] and (rule['family'] == 'A' or row['B']['positive'] is True)


def prepare():
    require(not (RUN / 'plan.json').exists(), 'Prepared experiment exists; do not overwrite')
    info_plan = information.p.verify_inputs()  # Sources only, no old predictions/labels.
    sample = load(v4.RUN / 'plan.json')['candidates']
    freeze = load(v4.RUN / 'output_freeze.json')['files']
    current = {(r['record']['case_id'], r['candidate']['organism_key']): r for r in v4.all_rows()}
    input_index = {(r['case_id'], r['organism_key']): r for r in info_plan['candidates']}
    rows, jobs = [], []
    for s in sample:
        row = current[(s['case_id'], s['organism_key'])]
        source = information.RUN / input_index[(s['case_id'], s['organism_key'])]['folder']
        folder = RUN / s['folder']
        core = load(source / 'core_packet.json'); view = a1_view(core)
        save(folder / 'a1_patient_view.json', view)
        write_new(folder / 'patient_provenance_not_model_input.json', (source / 'provenance_not_model_input.json').read_bytes())
        for name in ['pubmed.json', 'pubmed_raw.xml']:
            old_file = v4.RUN / s['folder'] / name
            if name == 'pubmed_raw.xml' and not old_file.exists():
                require(not load(v4.RUN / s['folder'] / 'pubmed.json')['articles'], 'Missing XML for nonempty literature')
                continue  # A zero-hit query has no efetch XML; don't invent one.
            require(digest(old_file.read_bytes()) == freeze[s['folder'] + '/' + name], 'Frozen literature changed')
            write_new(folder / name, old_file.read_bytes())
        articles = load(folder / 'pubmed.json')['articles']
        require(len(articles) <= 10 and len({a['pmid'] for a in articles}) == len(articles), 'Invalid frozen article list')
        save(folder / 'articles_with_sentences.json', [article_with_sentences(a) for a in articles])
        branch_data = row['candidate']['analytical_evidence']['per_test_signals']
        b = b_strength(branch_data); save(folder / 'B_not_model_input.json', b)
        rows.append({k: s[k] for k in ['case_id', 'group', 'organism', 'organism_key', 'baseline_tier', 'candidate_token', 'folder']})
        rows[-1].update({'article_count': len(articles), 'abstract_count': sum(bool(a['abstract']) for a in articles),
            'source_patient_result': row['record']['patient_result_file'],
            'source_patient_result_sha256': row['record']['patient_result_sha256'], 'B': b})
        jobs.extend({'case_id': s['case_id'], 'folder': s['folder'], 'pmid': a['pmid'], 'replicate': 0,
                     'has_abstract': bool(a['abstract'])} for a in articles)
    # Six first-article repeats, two cases per group; source-order, no label access.
    for group in ['KH33', 'KH5', 'A11']:
        for row in [r for r in rows if r['group'] == group and r['abstract_count']][:2]:
            main = next(j for j in jobs if j['case_id'] == row['case_id'] and j['has_abstract'])
            jobs.append({**main, 'replicate': 1})
    dependencies = [Path(__file__), Path(old.__file__), Path(v4.__file__), HERE / 'inspect_release.py', HERE / 'api_probe.py',
        Path(information.__file__), Path(information.p.__file__), LEGACY / 'rag_re_casefit/rules.py', SOURCE_PROMPT,
        information.RUN / 'input_freeze.json', information.RUN / 'plan.json', v4.RUN / 'output_freeze.json']
    plan = {'objective': 'Micro F1; precision decrease permitted; no patient-specific exceptions',
        'selection_reads_answers': False, 'independent_test': False, 'full49_candidate_run': False,
        'sample': 'Same fixed 18 candidate-case pairs as V4/V5; 16 Context +2 existing positive controls. Previously viewed development data.',
        'model': MODEL, 'reasoning': 'medium', 'max_output_tokens': 3500, 'max_workers': 4,
        'legacy_A1_truth_table_unchanged': True, 'patient_fact_view_changed_to_verified_sources': True,
        'candidate_level_LLM_recommendation': False, 'new_pubmed_queries': False,
        'primary_rule_ids': ['E0', 'B_only', 'A_strong_k2', 'AB_strong_k2'],
        'primary_k': 2, 'k_reason': 'Historical strong>=2 reference, fixed before this experiment, not selected from current labels',
        'secondary_match_reference_k': 7, 'all_rules': rule_configs(),
        'B_definition': 'B_v7_relative_rank3; SAME branch reads>0, QC evaluable/partial_evaluable, actual rank_in_retained_universe<=3. Not old B_STRICT, not a validated diagnostic threshold.',
        'B_rank_limit': B_RANK_LIMIT, 'accepted_qc': ACCEPTED_QC,
        'missing_abstract_action': 'insufficient, no API call; never pad to ten or infer missing abstract',
        'counterevidence_action': 'Report context, do not automatically veto this patient due to another published case.',
        'primary_counts_use_replicate': 0, 'repeat_selection': 'first usable article, first two cases per group in frozen source order',
        'dependencies': {str(p.relative_to(ROOT)).replace('\\','/'): digest(p.read_bytes()) for p in dependencies},
        'prompt_sha256': digest(PROMPT.encode()), 'schema_sha256': digest(old.encoded(SCHEMA)), 'candidates': rows, 'jobs': jobs}
    save(RUN / 'plan.json', plan)
    write_new(RUN / 'program_snapshot.py', Path(__file__).read_bytes())
    write_new(RUN / 'prompt_snapshot.txt', PROMPT.encode())
    files = {p.relative_to(RUN).as_posix(): digest(p.read_bytes()) for p in RUN.rglob('*') if p.is_file()}
    save(RUN / 'input_freeze.json', {'answers_read': False, 'files': files})
    verify_inputs()
    print(json.dumps({'candidates': len(rows), 'main_articles': sum(r['article_count'] for r in rows),
        'main_abstracts': sum(r['abstract_count'] for r in rows), 'repeat_calls': sum(j['replicate'] == 1 for j in jobs),
        'B_true': sum(r['B']['positive'] is True for r in rows), 'B_unknown': sum(r['B']['positive'] is None for r in rows)}, ensure_ascii=True), flush=True)


def verify_inputs():
    plan = load(RUN / 'plan.json')
    for name, sha in load(RUN / 'input_freeze.json')['files'].items():
        require(digest((RUN / name).read_bytes()) == sha, 'Frozen input changed: ' + name)
    for name, sha in plan['dependencies'].items():
        require(digest((ROOT / name).read_bytes()) == sha, 'Dependency changed: ' + name)
    information.p.verify_inputs()
    for r in plan['candidates']:
        require(digest((HERE / r['source_patient_result']).read_bytes()) == r['source_patient_result_sha256'], 'Candidate source changed')
    return plan


def request_for(row, view, article):
    return {'model': MODEL, 'reasoning': {'effort': 'medium'}, 'store': False, 'max_output_tokens': 3500,
        'instructions': PROMPT, 'input': json.dumps({'organism': row['organism'],
            'documented_clinical_names_not_scoring_aliases': CLINICAL_NAMES.get(row['organism'], []),
            'patient_case_card': view, 'article': article}, ensure_ascii=False),
        'text': {'format': {'type': 'json_schema', 'name': 'ober_article_fields_only', 'strict': True, 'schema': SCHEMA}}}


def job(item, row, env):
    folder = RUN / item['folder']; view = load(folder / 'a1_patient_view.json')
    article = next(a for a in load(folder / 'articles_with_sentences.json') if a['pmid'] == item['pmid'])
    target = folder / 'article_judgments' / f"{item['pmid']}_r{item['replicate']}"
    if not article['abstract']:
        save(target.with_suffix('.json'), no_abstract(view, article)); return
    request = request_for(row, view, article); sha = digest(old.encoded(request))
    save(target.parent / (target.name + '_request.json'), request)
    raw_file = target.parent / (target.name + '_raw_response.json')
    if raw_file.exists():
        record = load(raw_file); require(record['request_sha256'] == sha, 'Stale response cache')
        raw = record['response']
    else:
        base = env.get('OPENAI_BASE_URL', 'https://api.openai.com/v1').rstrip('/')
        require(base == 'https://api.openai.com/v1', 'Unapproved API endpoint')
        started = time.monotonic()
        raw = json.loads(old.http(base + '/responses', request, env['OPENAI_API_KEY'], timeout=240))
        save(raw_file, {'request_sha256': sha, 'elapsed_seconds': time.monotonic()-started, 'response': raw})
    require(raw.get('status') == 'completed', 'Incomplete API response retained; not a negative result and not auto-retried')
    text = ''.join(c.get('text','') for o in raw.get('output',[]) for c in o.get('content',[]) if c.get('type') == 'output_text')
    output = json.loads(text)
    checked = validate_extraction(output, view, article)
    save(target.with_suffix('.json'), {'status': 'completed', 'pmid': article['pmid'], 'replicate': item['replicate'],
        'request_sha256': sha, 'model_returned': raw.get('model'), 'response_id': raw.get('id'), 'usage': raw.get('usage'), **checked})
    print(json.dumps({'done': item['case_id'], 'pmid': item['pmid'], 'replicate': item['replicate'],
        'verdict': checked['verdict'], 'trace_valid': checked['traceability_valid']}), flush=True)


def derive_decisions(plan):
    features = []
    for row in plan['candidates']:
        folder = RUN / row['folder']
        articles = load(folder / 'articles_with_sentences.json')
        judgments = [load(folder / 'article_judgments' / (a['pmid'] + '_r0.json')) for a in articles]
        aggregate = aggregate_judgments(judgments)
        save(folder / 'A1_aggregate.json', aggregate)
        features.append({**row, 'A1': aggregate, 'invalid_trace_count': sum(not j['traceability_valid'] for j in judgments)})
    selections = [{'rule': rule, 'selected_tokens': [r['candidate_token'] for r in features if admit(r, rule)],
                   'selected': [{'case_id': r['case_id'], 'organism': r['organism']} for r in features if admit(r, rule)]}
                  for rule in plan['all_rules']]
    save(RUN / 'decisions_without_answers.json', {'features': features, 'selections': selections,
        'interpretation': 'Count gates are research decisions. A1 case fit and B analytical strength do not establish causal infection.'})


def verify_outputs():
    freeze = load(RUN / 'output_freeze.json')
    for name, sha in freeze['files'].items():
        require(digest((RUN / name).read_bytes()) == sha, 'Frozen output changed: ' + name)


def run(smoke=False):
    plan = verify_inputs()
    if (RUN / 'output_freeze.json').exists():
        verify_outputs(); print('Frozen run verified; no API or credentials.', flush=True); return
    env = environment(); require(bool(env.get('OPENAI_API_KEY')), 'No existing credential')
    rows = {r['case_id']: r for r in plan['candidates']}
    jobs = [j for j in plan['jobs'] if j['replicate'] == 0 and j['has_abstract']][:2] if smoke else plan['jobs']
    errors = []
    def work(item):
        try:
            job(item, rows[item['case_id']], env)
        except Exception as exc:
            errors.append({'case': item['case_id'], 'pmid': item['pmid'], 'replicate': item['replicate'], 'error_type': type(exc).__name__})
            print(json.dumps({'failed': errors[-1]}), flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=plan['max_workers']) as executor:
        list(executor.map(work, jobs))
    require(not errors, 'Technical failures, not negative decisions: ' + json.dumps(errors))
    if smoke:
        print('Two independent article calls completed; no reference access.', flush=True); return
    derive_decisions(plan)
    paths = [p for p in (RUN / 'cases').rglob('*') if p.is_file()] + [RUN / 'decisions_without_answers.json']
    save(RUN / 'output_freeze.json', {'answers_read': False,
        'files': {p.relative_to(RUN).as_posix(): digest(p.read_bytes()) for p in paths},
        'input_freeze_sha256': digest((RUN / 'input_freeze.json').read_bytes())})
    print('All article extractions AND all k decisions frozen before reference access.', flush=True)


def metrics(tp, fp, fn):
    return v4.metrics(tp, fp, fn)


def evaluate():
    plan = verify_inputs(); verify_outputs()
    # First reference access in this runner, after all article and decision freeze.
    from inspect_release import diagnostic_name
    reference = load(HERE / 'cohort49_review.json')
    truths = {r['case_id']: {diagnostic_name(n) for n in r['reference']} for r in reference['cases']}
    baseline = {r['case_id']: {diagnostic_name(n) for n in r['predicted']} for r in reference['cases']}
    decision = load(RUN / 'decisions_without_answers.json')
    features = [{**f, 'matches_reference': diagnostic_name(f['organism']) in truths[f['case_id']]} for f in decision['features']]
    by_token = {r['candidate_token']: r for r in features}
    scores = []
    for selected in decision['selections']:
        added = [by_token[t] for t in selected['selected_tokens']]
        predictions = {cid: set(names) for cid,names in baseline.items()}
        for row in added:
            predictions[row['case_id']].add(diagnostic_name(row['organism']))
        group_scores = {}
        for group in ['ALL49','KH33','KH5','A11']:
            ids = {r['case_id'] for r in reference['cases'] if group == 'ALL49' or r['analysis_group'] == group}
            m = metrics(sum(len(predictions[c] & truths[c]) for c in ids), sum(len(predictions[c]-truths[c]) for c in ids), sum(len(truths[c]-predictions[c]) for c in ids))
            group_added = [r for r in added if r['case_id'] in ids]
            tp = sum(r['matches_reference'] and diagnostic_name(r['organism']) not in baseline[r['case_id']] for r in group_added)
            fp = sum(not r['matches_reference'] for r in group_added)
            group_scores[group] = {'metrics': m, 'added_tp': tp, 'added_fp': fp,
                'delta_f1': m['f1'] - reference['groups'][group]['f1']}
        scores.append({'rule': selected['rule'], 'selected': selected['selected'], 'groups': group_scores})
    e0 = next(s for s in scores if s['rule']['id'] == 'E0')['groups']['ALL49']['metrics']
    require(all(e0[k] == reference['groups']['ALL49'][k] for k in ['tp','fp','fn']), 'Baseline contract mismatch')
    repeats = []
    usage = Counter()
    for j in plan['jobs']:
        folder = RUN / j['folder'] / 'article_judgments'
        output = load(folder / f"{j['pmid']}_r{j['replicate']}.json")
        if output['usage']:
            for k in ['input_tokens','output_tokens','total_tokens']:
                usage[k] += output['usage'].get(k, 0)
        if j['replicate'] == 1:
            first = load(folder / (j['pmid'] + '_r0.json'))
            repeats.append({'case_id': j['case_id'], 'pmid': j['pmid'], 'same_full_request': first['request_sha256'] == output['request_sha256'],
                'same_verdict': first['verdict'] == output['verdict'], 'first': first['verdict'], 'repeat': output['verdict'],
                'field_agreement': {k: first['extraction'][k] == output['extraction'][k] for k in FIELDS}})
    result = {'primary_rule_ids': plan['primary_rule_ids'], 'not_independent_test': True, 'not_full49_candidate_run': True,
        'baseline': reference['groups'], 'scores': scores, 'features': features, 'repeat_checks': repeats, 'usage': dict(usage),
        'scoring_contract': 'Existing retrospective diagnostic_name reconciliation including approved group aliases; not pure exact species; Candida species remain distinct.'}
    save(RUN / 'evaluation.json', result)
    write_report(plan, result)
    print(json.dumps({'primary': [s for s in scores if s['rule']['id'] in plan['primary_rule_ids']],
        'usage': dict(usage), 'repeat_flips': sum(not r['same_verdict'] for r in repeats),
        'field_trace_invalid_articles': sum(f['invalid_trace_count'] for f in features)}, ensure_ascii=True), flush=True)


def write_report(plan, result):
    lines = ['# OBER：逐篇A1＋固定補入規則', '',
        '18個固定開發病例－候選，16個Context＋2個既有陽性對照。不是test，也不是全49例候選重跑。',
        '原97筆陽性不動。下表為只加入本輪候選後的49例反事實micro分數；答案／名稱配對契約不改，不同Candida species不合併。',
        'LLM只抽取逐篇欄位與雙邊來源，不輸出／使用整體補入建議。舊A1 truth table未改，程式計strong/partial，再以固定條件補入。',
        'A1不看mNGS數值、culture/PCR結果、上游選菌理由或醫師答案。背景／影像沿用已核對核心來源，不再由LLM重寫病例。', '',
        '## B的定義與限制', '',
        '新版缺舊B的signal tier、reads tier、dominance、percentile，故舊B_STRICT不可直接重現。本輪B_v7只測相對排序：同一分支reads>0、QC為evaluable或partial_evaluable、rank_in_retained_universe<=3。',
        '這不是原始全菌宇宙的percentile，不用analytical_rank_priority冒充實際rank，不跨分支加總reads，不把多分支當獨立確認；RPM只記錄。不是經臨床驗證的感染門檻。', '',
        '## 預先指定主比較：A1 strong >=2', '',
        '| 規則 | 新增TP | 新增FP | Precision | Recall | F1 |', '|---|---:|---:|---:|---:|---:|']
    for s in result['scores']:
        if s['rule']['id'] in plan['primary_rule_ids']:
            g = s['groups']['ALL49']; m = g['metrics']
            lines.append(f'| {s["rule"]["id"]} | {g["added_tp"]} | {g["added_fp"]} | {m["precision"]:.2%} | {m["recall"]:.2%} | {m["f1"]:.2%} |')
    lines += ['', 'E0＝原清單；A＝原清單OR文獻門檻；B＝原清單OR B_v7；AB＝原清單OR(A門檻AND B_v7)。所有補入只處理Context，對照不刪原陽性。',
        'AB的補入一定是B的子集合，不能宣稱它比B單獨使用提高recall；文獻過濾是否值得要看保留TP／排除FP後的F1。', '',
        '## k完整曲線（不是挑最佳點當驗證）', '',
        '| 計數 | k | A新增TP/FP | A的P/R/F1 | AB新增TP/FP | AB的P/R/F1 |', '|---|---:|---|---|---|---|']
    index = {s['rule']['id']: s for s in result['scores']}
    for count in ['strong','match']:
        for k in K_VALUES:
            ga = index[f'A_{count}_k{k}']['groups']['ALL49']; gb = index[f'AB_{count}_k{k}']['groups']['ALL49']
            def fmt(g):
                m = g['metrics']; return f'{m["precision"]:.2%} / {m["recall"]:.2%} / {m["f1"]:.2%}'
            lines.append(f'| {count} | {k} | {ga["added_tp"]}/{ga["added_fp"]} | {fmt(ga)} | {gb["added_tp"]}/{gb["added_fp"]} | {fmt(gb)} |')
    lines += ['', 'strong＝舊A1 strong_match篇數；match＝strong+partial篇數。strong k=2是歷史參照而非本輪最佳化；match k=7只是既有參照／敏感度分析。',
        '沒有把舊aggregate的「high需至少5篇judgeable」額外塞入篇數門檻；本輪測的是舊rescue使用的strong/match計數。aggregate tier仍保留供查核，不當另一個隱藏閘門。', '',
        '## 逐候選與病原分布', '',
        '| 病例 | 病原 | 原層級 | 答案含此菌 | 文獻/有摘要 | strong | partial | 不足 | B_v7 | A2 | AB2 |',
        '|---|---|---|---|---|---:|---:|---:|---|---|---|']
    for f in result['features']:
        a = f['A1']
        lines.append(f'| {f["case_id"]} | {f["organism"]} | {f["baseline_tier"]} | {f["matches_reference"]} | {f["article_count"]}/{f["abstract_count"]} | {a["strong_match_count"]} | {a["partial_match_count"]} | {a["insufficient_count"]} | {f["B"]["positive"]} | {admit(f, index["A_strong_k2"]["rule"])} | {admit(f, index["AB_strong_k2"]["rule"])} |')
        folder = RUN / f['folder']; detail = [f'# {f["case_id"]} — {f["organism"]}', '',
            '逐篇資料抽取，不是LLM整體診斷。此人讀檔在輸出凍結及評分後產生。',
            f'固定答案含此候選：{f["matches_reference"]}；原層級：{f["baseline_tier"]}；B_v7：{f["B"]["positive"]}。',
            f'A1 strong={a["strong_match_count"]}、partial={a["partial_match_count"]}、不足={a["insufficient_count"]}。', '', '## 逐篇']
        for article in load(folder / 'articles_with_sentences.json'):
            j = load(folder / 'article_judgments' / (article['pmid']+'_r0.json'))
            detail += ['', f'### PMID {article["pmid"]} — {j["verdict"]}',
                f'https://pubmed.ncbi.nlm.nih.gov/{article["pmid"]}/', j['extraction']['rationale'],
                '欄位：' + json.dumps({key:j['extraction'][key] for key in FIELDS}, ensure_ascii=False),
                '來源檢查：' + str(j['traceability_valid'])]
        write_new(folder / 'finding_zh.md', ('\n'.join(detail)+'\n').encode())
    lines += ['', '## 重複判讀與用量', '',
        f'六個預先指定文章配對，同請求重抽的verdict翻轉：{sum(not r["same_verdict"] for r in result["repeat_checks"])} / {len(result["repeat_checks"])}。主分數一律使用第一次，不挑較好那次。',
        '```json', json.dumps(result['usage'], indent=2), '```', '',
        '## 解讀限制', '',
        '- 文獻匹配不是病人感染機率，十篇不是十次獨立診斷；相同病例或重疊研究仍可能跨論文重複，本輪只保證PMID不重複。',
        '- 無摘要／找不到文獻保留不足，不補滿十篇、不當病原不致病。文獻庫與檢索覆蓋會限制A。',
        '- 一篇定植或再活化文獻不自動否定這位病人；只保留反向情境，不加入一票否決規則。',
        '- 來源ID和短引文可定位不代表語意或醫學判讀已驗證；trace failure會成insufficient，不當臨床陰性。',
        '- 未知或跨事件影像不可當已證實當次phenotype；史料矛盾保持未知。',
        '- strong的額外匹配仍可能只來自同類檢體；舊定義保留作基準，不宣稱已解決病例區辨。',
        '- 只看開發資料最高F1會樂觀。完整k曲線不等於有獨立test；不依特定病人或菌種寫例外。',
        '- 此輪沒有C3、二次網路檢索、自主擴窗或最後LLM自由補入。',
        '- 名稱配對沿用既有核對契約，包含核准group aliases，並非pure exact species；Candida各species分開。', '',
        'API介面：[OpenAI Docs Structured outputs](https://developers.openai.com/api/docs/guides/structured-outputs)。結構正確不代表證據內容正確。']
    write_new(RUN / 'REPORT_zh.md', ('\n'.join(lines)+'\n').encode())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('command', choices=['prepare','verify','smoke','run','evaluate'])
    args = parser.parse_args()
    {'prepare': prepare, 'verify': verify_inputs, 'smoke': lambda: run(True), 'run': run, 'evaluate': evaluate}[args.command]()
