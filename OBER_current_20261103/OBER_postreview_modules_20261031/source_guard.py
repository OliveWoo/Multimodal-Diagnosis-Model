"""Answer-blind sensitivity guard for already-documented temporal input defects."""
import collections
from pathlib import Path
import modules as m


def main():
    m.verify('prediction');m.require(not (m.RUN/'source_guard_plan.json').exists(),'Guard already frozen')
    prior=m.ROOT/'OBER_case_lookup_20261031/runs/lookup_v1'
    for path,h in m.load(prior/'source_freeze.json')['files'].items():
        m.require(m.digest(Path(path).read_bytes())==h,'Changed previous source audit')
    audited={r['token']:r for r in m.load(prior/'source_plan.json')['candidates']}
    rows=[];counts=collections.Counter();files=[]
    for row in m.load(m.RUN/'predictions_without_answers.json'):
        a=audited[row['token']];folder=prior/a['folder']
        original=m.load(folder/'control_original.json');safe=m.load(folder/'control_time_checked.json')
        frozen=m.load(m.PARENT/'judgements'/row['case_id']/row['token']/'T/input.json')
        m.require(original==frozen,'Temporal audit refers to different original input')
        old={f['id'] for f in original['patient_facts']};new={f['id'] for f in safe['patient_facts']}
        missing=sorted(old-new);m.require(bool(missing)==a['control_time_repaired'],'Temporal audit disagreement')
        decisions={s:bool(v and not missing) for s,v in row['decisions'].items()}
        rows.append({'token':row['token'],'case_id':row['case_id'],'deferred_due_to_known_later_input':bool(missing),
            'known_later_fact_ids':missing,'decisions':decisions})
        counts['candidate_inputs_with_known_later_records']+=bool(missing)
        counts['T_broad_proposals_with_known_later_input']+=bool(missing) and row['features']['T_broad']
        files += [folder/'control_original.json',folder/'control_time_checked.json']
    m.save(m.RUN/'source_guard_predictions_without_answers.json',rows)
    m.save(m.RUN/'source_guard_plan.json',{'no_gold_loaded':True,'no_new_API_calls':True,
        'same_guard_all24_settings':True,'no_patient_or_reference_removed':True,'counts':dict(counts),
        'previous_audit_freeze_sha256':m.digest((prior/'source_freeze.json').read_bytes())})
    m.freeze('source_guard',files+[m.RUN/'source_guard_predictions_without_answers.json',m.RUN/'source_guard_plan.json',
        Path(__file__),m.HERE/'SOURCE_AMENDMENT_01_zh.md'])
    print(dict(counts),flush=True)


if __name__=='__main__':main()
