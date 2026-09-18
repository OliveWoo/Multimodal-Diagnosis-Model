"""Run existing unittest tests and plain test functions with tmp_path fixtures."""
from pathlib import Path
import sys,unittest,inspect,tempfile,importlib

root=Path(sys.argv[1]).resolve()
testdir=root/(sys.argv[2] if len(sys.argv)>2 else 'tests')
exclude=set(sys.argv[3:])
sys.path.insert(0,str(root));sys.path.insert(0,str(testdir))
suite=unittest.TestSuite()
for p in sorted(testdir.glob('test*.py')):
    if p.name in exclude:continue
    module=importlib.import_module(p.stem)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromModule(module))
    for name,func in inspect.getmembers(module,inspect.isfunction):
        if name.startswith('test_') and func.__module__==module.__name__:
            params=list(inspect.signature(func).parameters)
            if params not in ([],['tmp_path']):raise RuntimeError(f'Unsupported fixture in {p}:{name}: {params}')
            def wrapped(f=func,params=params):
                if params:
                    with tempfile.TemporaryDirectory() as t:f(Path(t))
                else:f()
            suite.addTest(unittest.FunctionTestCase(wrapped,description=p.name+':'+name))
result=unittest.TextTestRunner(verbosity=2).run(suite)
raise SystemExit(0 if result.wasSuccessful() else 1)

