"""Isolated A1-on-S artifacts, reusing unchanged legacy grading."""
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'OBER_postreview_modules_20261031'))
import modules as source
import literature_rule_pilot as legacy

c=source.old
RUN=HERE/'runs/a1s_v1'
load,save,require,digest,write_new=c.load,c.save,c.require,c.digest,c.write_new
MODEL=c.MODEL
PRIMARY='T_ONLY_strong_k2_BOTH'


def freeze(stage,paths):
    save(RUN/(stage+'_freeze.json'),{'files':{str(p.resolve()):digest(p.read_bytes()) for p in paths},
        'posthoc_development_not_independent_test':True})


def verify(stage):
    for path,h in load(RUN/(stage+'_freeze.json'))['files'].items():
        require(digest(Path(path).read_bytes())==h,'Changed frozen A1/S artifact '+path)


def configs():
    return [{'id':f'{scope.upper()}_{count}_k{k}_{signal.upper()}',
             'scope':scope,'count':count,'k':k,'signal':signal}
        for scope in ['T_only','BOTH_counts'] for count in ['strong','match']
        for k in range(11) for signal in ['none','either','both']]


def admit(row,counts,rule):
    if not row['S_pass']:return False
    if rule['signal']!='none' and not source.native.passes(row['native_flags'],rule['signal']):return False
    if rule['k']==0:return True
    require(counts is not None,'Missing A1 is unavailable, not a count of zero')
    def value(arm):
        x=counts[arm];return x['strong_match_count']+(x['partial_match_count'] if rule['count']=='match' else 0)
    score=value('T') if rule['scope']=='T_only' else min(value('G'),value('T'))
    return score>=rule['k']
