"""Entry point with a pre-run correction to reporting-only scope-key mapping."""
import argparse
from pathlib import Path
import repeat_driver as d
import compare_runs as comparison


def profile(record):
    if record['status']=='unresolved_qc':return ['unresolved_qc']
    # normalize() stores 'scope'; human-readable signature columns call it
    # 'organism_scope'. Do not silently compare None for all resolved identities.
    return [record['status'],record.get('scope',record.get('organism_scope'))]+[
        record.get(k) for k in comparison.SIGNATURE_FIELDS[1:]]


def bind():
    path=d.BUNDLE/'comparison_mapping_freeze.json'
    if path.exists():d.verify_file_manifest(path)
    else:d.freeze_bundle('comparison_mapping',[Path(__file__),d.HERE/'REPORT_MAPPING_zh.md',d.HERE/'test_profile_mapping.py'])
    comparison.profile=profile


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['pipeline','evaluate'])
    bind();{'pipeline':d.pipeline,'evaluate':d.evaluate}[p.parse_args().stage]()
