"""Explicit new run for incomplete replayed tool responses; not a silent retry.

Keeps identical clinical facts/judge. Each lookup planner reads an explicit fresh
snapshot of the currently visible facts rather than replaying opaque API history.
Reuses r2 completed responses ONLY when the full request hash is identical.
"""
import argparse
import json
from pathlib import Path

import information_pilot as p

PREVIOUS = p.RUN
RUN = p.HERE / 'runs/pilot_v5_information_20261031_r3'
p.RUN = RUN
original_api_call = p.api_call


def prepare():
    p.require(not (RUN / 'plan.json').exists(), 'Recovery plan already exists')
    previous = p.load(PREVIOUS / 'plan.json')
    freeze = p.load(PREVIOUS / 'input_freeze.json')['files']
    for name, sha in freeze.items():
        p.require(p.digest((PREVIOUS / name).read_bytes()) == sha, 'Original frozen input changed')
        if name != 'plan.json':
            p.write_new(RUN / name, (PREVIOUS / name).read_bytes())
    previous['dependencies'][str(Path(__file__).relative_to(p.ROOT)).replace('\\', '/')] = p.digest(Path(__file__).read_bytes())
    previous['transport_recovery'] = {
        'reason': 'r2 returned incomplete/max_messages on replayed second-round tool responses; these are technical failures, not negative diagnoses.',
        'change': 'Fresh visible-fact snapshot for each lookup planner; same clinical source universe, max rounds, tool, judge and policies.',
        'cache': 'Reuse completed prior response only with identical full request SHA256.',
        'billing_caveat': 'Two in-flight requests at explicit r2 process stop may lack locally recorded response/usage; saved token sums are not a complete billing ledger.'}
    p.save(RUN / 'plan.json', previous)
    p.write_new(RUN / 'recovery_program_snapshot.py', Path(__file__).read_bytes())
    p.save(RUN / 'input_freeze.json', {'answers_read': False,
        'files': {q.relative_to(RUN).as_posix(): p.digest(q.read_bytes()) for q in RUN.rglob('*') if q.is_file()}})
    p.verify_inputs()
    print('Recovery input freeze verified; no label access or new API calls.', flush=True)


def cached_api_call(folder, name, request, env):
    local = folder / (name + '_raw_response.json')
    prior = PREVIOUS / folder.relative_to(RUN) / local.name
    if not local.exists() and prior.exists():
        record = p.load(prior)
        if record['request_sha256'] == p.digest(p.old.encoded(request)) and record['response'].get('status') == 'completed':
            p.write_new(local, prior.read_bytes())
            p.save(folder / (name + '_cache_provenance.json'), {
                'reused_from': str(prior.relative_to(p.ROOT)).replace('\\', '/'),
                'full_request_sha256': record['request_sha256'], 'new_api_charge': False})
    return original_api_call(folder, name, request, env)


def planner_history(full, visible, initial, log):
    data = {'case': p.packet_for(full, visible), 'initial_assessment': initial['assessment']}
    if log:
        data['prior_local_lookups'] = [{'round': q['round'], 'question': q['executed_arguments']['question'],
            'returned_fact_ids': [f['id'] for f in q['result']['facts']],
            'remaining_matching_facts': q['result']['remaining_matching_facts']} for q in log]
    return [{'role': 'user', 'content': json.dumps(data, ensure_ascii=False)}]


def run_case(row, env):
    folder = RUN / row['folder']; full = p.load(folder / 'source_catalog.json')
    visible = set(row['core_fact_ids']); core = p.packet_for(full, visible)
    initial = p.judge(folder, 'I0', core, env)
    log = []
    for round_no in range(p.MAX_ROUNDS):
        if len(visible) == len(full['facts']):
            break
        request = p.base_request(p.PLANNER, planner_history(full, visible, initial, log))
        request.update({'tools': [p.TOOL], 'tool_choice': {'type': 'function', 'name': p.TOOL['name']},
                        'parallel_tool_calls': False, 'include': ['reasoning.encrypted_content']})
        raw = p.api_call(folder, f'I1_lookup_{round_no+1}', request, env)
        p.require(raw.get('status') == 'completed', 'Incomplete tool response; not a negative decision')
        calls = [item for item in raw['output'] if item['type'] == 'function_call']
        p.require(len(calls) == 1 and calls[0]['name'] == p.TOOL['name'], 'Unexpected tool call')
        proposed = json.loads(calls[0]['arguments']); args = proposed
        result = p.lookup(full, visible, args); forced = False
        if args['scope'] == 'none' and round_no == 0 and initial['assessment']['recommendation'] != 'recommend_add':
            args = {'question': 'Safety check of unread source evidence before an absence-based non-addition.',
                    'scope': 'all', 'terms': [], 'temporal_roles': []}
            result = p.lookup(full, visible, args); forced = True
        visible.update(f['id'] for f in result['facts'])
        log.append({'round': round_no+1, 'call_id': calls[0]['call_id'], 'model_arguments': proposed,
            'executed_arguments': args, 'forced_safety_check': forced, 'result': result})
        if args['scope'] == 'none':
            break
    p.save(folder / 'I1_lookup_log.json', log)
    if visible == set(row['core_fact_ids']):
        p.save(folder / 'I1_visible_packet.json', core)
        p.save(folder / 'I1_result.json', {**initial, 'reused_I0_identical_information': True})
    else:
        p.judge(folder, 'I1', p.packet_for(full, visible), env)
    p.judge(folder, 'I2', p.packet_for(full, {f['id'] for f in full['facts']}), env)


def smoke():
    plan = p.verify_inputs(); env = p.environment()
    p.require(bool(env.get('OPENAI_API_KEY')), 'Missing API credential')
    for row in plan['candidates'][:2]:
        run_case(row, env)
    print('Two-case recovery smoke complete. No reference access; full outputs not yet frozen.', flush=True)


p.api_call = cached_api_call
p.run_case = run_case


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('command', choices=['prepare', 'smoke', 'run', 'verify', 'evaluate'])
    args = parser.parse_args()
    {'prepare': prepare, 'smoke': smoke, 'run': p.run, 'verify': p.verify_inputs, 'evaluate': p.evaluate}[args.command]()
