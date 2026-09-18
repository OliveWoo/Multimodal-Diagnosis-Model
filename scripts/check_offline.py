"""Run the portable, synthetic test suites without network or API credentials.

Usage: python scripts/check_offline.py [--node PATH] [--output DIRECTORY]
Private-cohort integration tests remain in source but are excluded here.
"""
from pathlib import Path
import argparse, json, os, re, shutil, subprocess, sys, tempfile

ROOT = Path(__file__).resolve().parents[1]
SUITES = [
    ('frontend', 'tests', []),
    ('scripts', 'tests', []),
    ('upstream', 'tests', []),
    ('RAG_re', 'tests', []),
    ('RAG_re_clinical', 'tests', []),
    ('RAG_re_casefit', 'tests', []),
    ('RAG_re_clinical_pruning', 'tests', []),
    ('RAG_re_clinical_pruning_abc', 'tests', ['test_locked_integration.py']),
    ('OBER_patient_delivery', 'tests', []),
    ('mngs_candidate_evidence_pipeline_20260902/mngs_candidate_evidence_pipeline', 'tests', []),
    ('LLM_test/direct_raw_benchmark', '.', ['test_gold_corrections_20260912.py', 'test_direct_raw_runner.py', 'test_input_identity.py']),
]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--node', default=shutil.which('node'))
    parser.add_argument('--output', type=Path, default=ROOT/'local_outputs/offline_checks')
    args = parser.parse_args()
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(PYTHONPATH=str(ROOT/'scripts/offline_guard'), PYTHONDONTWRITEBYTECODE='1', PYTHONIOENCODING='utf-8', PYTHONUTF8='1')
    results = []
    with tempfile.TemporaryDirectory(prefix='research-offline-') as temporary:
        env['TMP'] = env['TEMP'] = temporary
        for relative, testdir, exclude in SUITES:
            proc = subprocess.run([sys.executable, '-B', str(ROOT/'scripts/run_test_suite.py'), str(ROOT/relative), testdir, *exclude], cwd=ROOT/relative, env=env, capture_output=True, text=True, encoding='utf-8', timeout=180)
            logname = relative.replace('/', '__')+'.log'
            (output/logname).write_text(proc.stdout+'\n'+proc.stderr, encoding='utf-8')
            count = re.search(r'Ran (\d+) tests?', proc.stderr)
            skipped = re.search(r'OK \(skipped=(\d+)\)', proc.stderr)
            row = dict(suite=relative, exit_code=proc.returncode, tests=int(count.group(1)) if count else None, skipped=int(skipped.group(1)) if skipped else 0, excluded_private_tests=exclude, log=logname)
            results.append(row); print(json.dumps(row), flush=True)
        if args.node:
            proc = subprocess.run([args.node, '--test', 'tests/delivery.test.mjs'], cwd=ROOT/'OBER_patient_delivery', env=env, capture_output=True, text=True, encoding='utf-8', timeout=180)
            (output/'delivery_node.log').write_text(proc.stdout+'\n'+proc.stderr, encoding='utf-8')
            count = re.search(r'(?:#|ℹ) tests (\d+)', proc.stdout)
            results.append(dict(suite='delivery_node', exit_code=proc.returncode, tests=int(count.group(1)) if count else None, log='delivery_node.log'))
        else:
            results.append(dict(suite='delivery_node', exit_code=None, tests=None, status='NOT RUN: node executable unavailable'))
    (output/'summary.json').write_text(json.dumps(results, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    return 0 if all(r['exit_code']==0 for r in results) else 1

if __name__=='__main__':
    raise SystemExit(main())
