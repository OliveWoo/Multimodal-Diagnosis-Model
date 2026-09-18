"""Verify copied research files against SOURCE_MANIFEST.csv; no imports/API calls."""
from pathlib import Path
import csv,hashlib,json

def main():
    root=Path(__file__).resolve().parents[1]
    with (root/'SOURCE_MANIFEST.csv').open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    failures=[]
    for row in rows:
        path=(root/row['destination']).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            failures.append({'path':row['destination'],'reason':'missing or outside repository'});continue
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:
            failures.append({'path':row['destination'],'reason':'content differs from recorded source copy'})
    print(json.dumps({'source_files':len(rows),'verified':len(rows)-len(failures),'failures':failures},ensure_ascii=False,indent=2))
    return 1 if failures else 0

if __name__=='__main__':raise SystemExit(main())
