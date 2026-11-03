"""Read-only diagnostics of frozen A1 and raw B ingredients, no model/fix/eval run.

Writes only a NEW diagnostic artifact. Does not reconstruct a production B score,
select pathogens, tune thresholds, or alter any frozen input/output.
"""
import json
import math
import re
import unicodedata
from collections import Counter
from pathlib import Path

import literature_rule_pilot as pilot

ROOT = pilot.ROOT
OUT = pilot.HERE / 'audits/A1_B_causes_20261031/linked_event_verified'


def key(value):
    return re.sub('[^a-z0-9]', '', unicodedata.normalize('NFKC', str(value)).casefold())


def number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def code_priority(tags):
    codes = {re.sub(r'\s+', '', str(x)).upper().replace('=RK00','=RK0').replace('=RK08','=RK8') for x in tags}
    ntc, rk0, rk8 = 'CODE_NTC=00' in codes, 'CODE_RK_NTC=RK0' in codes, 'CODE_RK_NTC=RK8' in codes
    return 1 if ntc and rk0 else 2 if ntc and rk8 else 3 if rk0 else 4 if rk8 else 5 if ntc else None


def read_tier(value):
    if not number(value) or value <= 0: return None
    return 'R4_very_high' if value >= 10000 else 'R3_high' if value >= 1000 else 'R2_medium' if value >= 100 else 'R1_low' if value >= 10 else 'R0_trace'


def component_under_universe(rows, name, reads):
    if not rows: return None
    sorted_rows = sorted(rows, key=lambda r: (-r['reads'], key(r['name'])))
    top = sorted_rows[0]; second = sorted_rows[1]['reads'] if len(sorted_rows)>1 else 0
    ratio = top['reads']/max(second,1)
    dominance = ('D3_dominant' if ratio>=10 else 'D2_moderate' if ratio>=3 else 'D1_low') if key(top['name'])==key(name) else 'D0_not_top'
    return {'denominator': len(rows), 'reads_percentile_candidate': round(sum(r['reads']<=reads for r in rows)/len(rows),4),
            'dominance_candidate': dominance, 'top_name': top['name'], 'top_reads': top['reads'], 'second_reads': second,
            'not_authoritative_old_B': True}


def audit_a1(plan):
    rows, strong = [], []
    for candidate in plan['candidates']:
        folder = pilot.RUN / candidate['folder']
        aggregate = pilot.load(folder/'A1_aggregate.json')
        article_rows = []
        for article in pilot.load(folder/'articles_with_sentences.json'):
            result = pilot.load(folder/'article_judgments'/(article['pmid']+'_r0.json'))
            extraction = result['extraction']
            status = {k: extraction[k] for k in pilot.FIELDS}
            item = {'pmid': article['pmid'], 'verdict': result['verdict'], 'fields': status,
                'extra_matches': [k for k in ['specimen_alignment','host_alignment','phenotype_alignment','temporal_alignment'] if extraction[k]=='match'],
                'site_uses_protocol_target_only': extraction['field_evidence']['site_alignment']['patient_fact_ids']==['P_SITE'],
                'syndrome_uses_protocol_target_only': extraction['field_evidence']['syndrome_alignment']['patient_fact_ids']==['P_SYNDROME']}
            article_rows.append(item)
            if result['verdict']=='strong_match':
                strong.append({'case_id':candidate['case_id'],'organism':candidate['organism'], **item,
                    'rationale':extraction['rationale'],'field_evidence':extraction['field_evidence']})
        rows.append({'case_id':candidate['case_id'],'organism':candidate['organism'],
            'strong':aggregate['strong_match_count'],'partial':aggregate['partial_match_count'],'articles':article_rows})
    return {'candidates':rows,'strong_judgments':strong,
            'strong_extra_match_patterns': dict(Counter('+'.join(r['extra_matches']) for r in strong)),
            'strong_with_host_match':sum(r['fields']['host_alignment']=='match' for r in strong),
            'strong_with_temporal_match':sum(r['fields']['temporal_alignment']=='match' for r in strong)}


def audit_b():
    manifests = pilot.load(pilot.HERE/'cohort49_manifest.json')['cases']
    raw_cache, raw_tests, specimens, scopes, test_provenance = {}, {}, {}, {}, {}
    for rec in manifests:
        source = pilot.information.p.raw_file(rec); raw = pilot.load(source)
        data = raw['files']['input/mngs_raw.json']; raw_cache[rec['case_id']] = data
        target = rec['adjacent_specimen_code']
        matches = [c for c in data['cases'] if key(c.get('specimen_code'))==key(target) or key(c.get('case_review_id'))==key(target)]
        pilot.require(len(matches)==1,'Ambiguous raw specimen identity: '+rec['case_id'])
        specimens[rec['case_id']] = matches[0]
        # Do not infer event identity from patient number alone. Index specimen
        # and explicitly linked release events are separately bound to raw IDs.
        allowed = {key(target): None}
        for event in rec.get('release_linked_events') or []:
            event_key = key(event['case_review_id'])
            seqs = set(event['seq_ids'])
            if event_key in allowed and allowed[event_key] is not None:
                allowed[event_key].update(seqs)
            else:
                allowed[event_key] = seqs
        for event_key, permitted_seqs in allowed.items():
            event_cases = [c for c in data['cases'] if event_key in {key(c.get('specimen_code')),key(c.get('case_review_id'))}]
            pilot.require(len(event_cases)==1, 'Ambiguous linked raw event: '+rec['case_id']+'/'+event_key)
            event_case = event_cases[0]
            tests = [t for t in event_case['tests'] if permitted_seqs is None or t['seq_id'] in permitted_seqs]
            if permitted_seqs is not None:
                pilot.require({t['seq_id'] for t in tests}==permitted_seqs, 'Linked raw seq IDs incomplete: '+rec['case_id'])
            for test in tests:
                tkey = (rec['case_id'],test['seq_id'])
                if tkey in raw_tests:
                    pilot.require(raw_tests[tkey]==test and
                        test_provenance[tkey]['raw_case_review_id']==event_case.get('case_review_id'),
                        'Conflicting raw test identity: '+str(tkey))
                    continue
                raw_tests[tkey] = test
                test_provenance[tkey] = {'raw_specimen_code':event_case.get('specimen_code'),
                    'raw_case_review_id':event_case.get('case_review_id'),
                    'is_index_specimen':event_case is matches[0],
                    'binding':'index_specimen' if event_case is matches[0] else 'explicit_release_linked_event'}
        for organism in data['organisms']:
            for measurement in organism['measurements']:
                if (rec['case_id'],measurement.get('seq_id')) not in raw_tests: continue
                if not number(measurement.get('reads')) or measurement['reads']<=0 or measurement.get('raw_available') is False: continue
                tkey = (rec['case_id'],measurement['seq_id'])
                scopes.setdefault(tkey,[]).append({'name':organism['name'],'reads':measurement['reads'],
                    'selected':measurement.get('selected'),'tags':measurement.get('selection_tags') or []})
    results=[]; counts=Counter()
    for current in pilot.v4.all_rows():
        rec,candidate=current['record'],current['candidate']; cid=rec['case_id']
        raw=raw_cache[cid]; by_name=[o for o in raw['organisms'] if key(o['name'])==key(candidate['organism_name'])]
        for signal in candidate['analytical_evidence'].get('per_test_signals') or []:
            counts['reported_branches']+=1; seq=signal['seq_id']; test=raw_tests.get((cid,seq))
            matches=[m for o in by_name for m in o['measurements'] if m.get('seq_id')==seq]
            row={'case_id':cid,'organism':candidate['organism_name'],'seq_id':seq,
                'reported_reads':signal.get('reads'),'reported_ordinal_rank':signal.get('rank_in_retained_universe'),
                'reported_code_priority':signal.get('analytical_rank_priority'),'reported_rpm':signal.get('rpm_total')}
            if len(matches)!=1 or test is None:
                row['status']='identity_or_name_mapping_unresolved';row['raw_name_match_count']=len(by_name)
                row['raw_measurement_match_count']=len(matches);counts[row['status']]+=1;results.append(row);continue
            m=matches[0];counts['raw_branch_identity_bound']+=1
            row['raw_event_binding']=test_provenance[(cid,seq)]
            counts['index_specimen_branches' if row['raw_event_binding']['is_index_specimen'] else 'explicit_linked_other_specimen_branches']+=1
            row.update({'status':'bound','raw_reads':m['reads'],'raw_selection_tags':m.get('selection_tags') or [],
                'raw_selected':m.get('selected'),'reads_match':signal.get('reads')==m['reads'],
                'reads_tier_candidate':read_tier(m['reads']), 'code_priority_candidate':code_priority(m.get('selection_tags') or [])})
            counts['reads_match']+=row['reads_match']
            row['code_priority_match']=row['code_priority_candidate']==signal.get('analytical_rank_priority')
            counts['code_priority_match']+=row['code_priority_match']
            denom=(test.get('parameters') or {}).get('input_reads')
            row['input_reads_denominator']=denom
            if number(denom) and denom>0 and number(m['reads']):
                row['rpm_recomputed']=m['reads']/denom*1e6
                if number(signal.get('rpm_total')):
                    row['rpm_matches']=math.isclose(row['rpm_recomputed'],signal['rpm_total'],rel_tol=1e-10,abs_tol=1e-10)
                    counts['reported_numeric_rpm']+=1;counts['rpm_matches']+=row['rpm_matches']
                else:
                    counts['raw_rpm_possible_but_report_null']+=1
            elif number(signal.get('rpm_total')):counts['reported_rpm_without_raw_denominator']+=1
            else:counts['rpm_missing_both']+=1
            universe=scopes.get((cid,seq),[]);selected=[r for r in universe if r['selected'] is True]
            row['all_available_positive_universe']=component_under_universe(universe,candidate['organism_name'],m['reads'])
            row['selected_positive_universe']=component_under_universe(selected,candidate['organism_name'],m['reads'])
            if row['all_available_positive_universe'] and row['selected_positive_universe']:
                counts['universe_comparable']+=1
                counts['percentile_depends_on_universe']+=row['all_available_positive_universe']['reads_percentile_candidate']!=row['selected_positive_universe']['reads_percentile_candidate']
                counts['dominance_depends_on_universe']+=row['all_available_positive_universe']['dominance_candidate']!=row['selected_positive_universe']['dominance_candidate']
            results.append(row)
    return {'counts':dict(counts),'branches':results,
        'warning':'Ingredient feasibility only. Two alternative universes are diagnostic illustrations, NOT silently selected production denominators or an old-B reconstruction. No new admission rule/F1 evaluated.'}


def main():
    plan=pilot.verify_inputs();pilot.verify_outputs()
    a,b=audit_a1(plan),audit_b()
    sources=[pilot.HERE/'cohort49_manifest.json',pilot.RUN/'input_freeze.json',pilot.RUN/'output_freeze.json',
        ROOT/'RAG_re_casefit/rag_re_casefit/rules.py',ROOT/'RAG_re_clinical/rag_re_clinical/clinical_rules.py',
        ROOT/'research_repo_audit_20260917/cohort_revision_20260918/upstream/tools/deterministic_mngs_max_scorer_v19_recovered.py',
        ROOT/'research_repo_audit_20260917/clean_repo/upstream/tools/build_ranked_mngs_from_grouped_folder.py',
        ROOT/'research_repo_audit_20260917/clean_repo/upstream/tools/rank_mngs_candidate_microbes.py',
        ROOT/'research_repo_audit_20260917/clean_repo/upstream/tools/mngs_big_agent.py']
    result={'diagnostic_only':True,'new_api_calls':0,'production_rules_changed':False,'new_f1_evaluation':False,
        'A1':a,'B_ingredients':b,'source_hashes':{str(p.relative_to(ROOT)):pilot.digest(p.read_bytes()) for p in sources}}
    pilot.save(OUT/'diagnostic.json',result)
    print(json.dumps({'A1_strong_patterns':a['strong_extra_match_patterns'],'A1_strong_host_match':a['strong_with_host_match'],
        'A1_strong_temporal_match':a['strong_with_temporal_match'],'B_ingredient_counts':b['counts'],
        'unresolved_examples':[r for r in b['branches'] if r['status']!='bound'][:8]},ensure_ascii=True))


if __name__=='__main__':main()
