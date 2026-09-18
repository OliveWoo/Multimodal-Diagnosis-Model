"""Compare reconstructed source to its supplied private CPython 3.13 cache.

The bytecode is parsed only after SHA256 validation; it is never executed.
Filenames and line/debug positions are excluded. Bytecode, constants (including
nested code), symbol/local/closure tables, flags, signatures, stack size and
exception tables must match. This does not prove historical external inputs.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, marshal, sys, types

ROOT = Path(__file__).resolve().parents[1]


def code_signature(value):
    if isinstance(value, types.CodeType):
        return (value.co_code, tuple(code_signature(x) for x in value.co_consts),
                value.co_names, value.co_varnames, value.co_freevars, value.co_cellvars,
                value.co_flags, value.co_argcount, value.co_kwonlyargcount,
                value.co_posonlyargcount, value.co_exceptiontable, value.co_stacksize,
                value.co_name, value.co_qualname)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--component', required=True, choices=['scorer', 'summary'])
    parser.add_argument('--cache', required=True, type=Path)
    args = parser.parse_args()
    proof = json.loads((ROOT/'docs/UPSTREAM_RECOVERY_PROVENANCE.json').read_text(encoding='utf-8'))
    row = proof['components'][args.component]
    data = args.cache.read_bytes()
    if hashlib.sha256(data).hexdigest() != row['cache_sha256']:
        raise SystemExit('Cache SHA256 differs from the reviewed evidence')
    if sys.version_info[:2] != (3, 13) or data[:4] != importlib.util.MAGIC_NUMBER:
        raise SystemExit('Use CPython 3.13 with matching bytecode magic')
    source = (ROOT/row['destination']).read_text(encoding='utf-8-sig')
    old = marshal.loads(data[16:])
    new = compile(source, '<reconstructed-source>', 'exec')
    equal = code_signature(old) == code_signature(new)
    print(json.dumps({'component': args.component, 'module_code_equal': equal,
                      'functions': sum(isinstance(c, types.CodeType) for c in old.co_consts)}))
    return 0 if equal else 1


if __name__ == '__main__':
    raise SystemExit(main())
