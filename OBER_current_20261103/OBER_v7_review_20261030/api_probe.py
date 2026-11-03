"""Read existing project credentials without displaying or persisting secrets."""
import json
import os
import re
import urllib.error
import urllib.request
from pathlib import Path

ENV_FILE = Path(__file__).resolve().parents[1] / 'RAG/RAG/env.env'


def environment():
    values = dict(os.environ)
    if ENV_FILE.is_file():
        for raw in ENV_FILE.read_text(encoding='utf-8-sig').splitlines():
            line = raw.strip()
            if not line or line.startswith('#') or '=' not in line:
                continue
            key, value = line.split('=', 1)
            key = re.sub(r'^(?:export\s+|\$env:)', '', key.strip())
            value = value.strip().strip('"').strip("'")
            if re.fullmatch(r'[A-Z][A-Z0-9_]*', key) and key not in values:
                values[key] = value
    return values


def main():
    env = environment()
    key = env.get('OPENAI_API_KEY')
    print(json.dumps({'key_available': bool(key), 'credential_source': 'existing_project_env' if key else None}))
    if not key:
        return
    base = env.get('OPENAI_BASE_URL', 'https://api.openai.com/v1').rstrip('/')
    request = urllib.request.Request(base + '/models', headers={'Authorization': 'Bearer ' + key})
    try:
        with urllib.request.urlopen(request, timeout=25) as response:
            result = json.load(response)
        ids = [m['id'] for m in result.get('data', [])]
        print(json.dumps({'models_readable': True, 'relevant_models': sorted(m for m in ids if any(s in m for s in ['luna', 'astra', 'sol', 'terra']))}))
    except urllib.error.HTTPError as exc:
        print(json.dumps({'models_readable': False, 'http_status': exc.code}))
    except Exception as exc:
        print(json.dumps({'models_readable': False, 'error_type': type(exc).__name__}))


if __name__ == '__main__':
    main()
