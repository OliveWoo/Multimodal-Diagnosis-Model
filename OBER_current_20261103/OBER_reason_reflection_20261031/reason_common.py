"""Shared immutable artifacts, official API transport and bounded dispatcher."""
import concurrent.futures
import json
import sys
import time
import urllib.error
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
OLD=ROOT/'OBER_literature_reflection_20261031'
sys.path.insert(0,str(OLD))
import reflection as old

RUN=HERE/'runs/reason_v1'
SOURCE=OLD/'runs/reflection_linked_v2'
load,save,require,digest,write_new=old.load,old.save,old.require,old.digest,old.write_new
obj,enum,string=old.obj,old.enum,old.string
MODEL=old.MODEL
WORKERS=32


def ids(values):
    values=list(dict.fromkeys(values))
    return {'type':'array','items':enum(values) if values else string(),**({'maxItems':0} if not values else {})}


def freeze(stage,paths):
    save(RUN/(stage+'_freeze.json'),{'files':{str(p.resolve()):digest(p.read_bytes()) for p in paths},
        'not_independent_test':True})


def verify(stage):
    data=load(RUN/(stage+'_freeze.json'))
    for name,h in data['files'].items():require(digest(Path(name).read_bytes())==h,'Changed frozen file '+name)
    return data


def deps(*names):
    return [HERE/'reason_common.py',HERE/'PROTOCOL_zh.md']+[HERE/n for n in names]


def api(folder,tag,request,env):
    sha=digest(old.d.clinical.sources.old.encoded(request));raw=folder/(tag+'_raw_response.json')
    save(folder/(tag+'_request.json'),request)
    if raw.exists():
        record=load(raw);require(record['request_sha256']==sha,'Stale API cache');return record['response']
    base=env.get('OPENAI_BASE_URL','https://api.openai.com/v1').rstrip('/')
    require(base=='https://api.openai.com/v1','Only official endpoint permitted')
    for attempt in range(3):
        try:
            start=time.monotonic()
            data=json.loads(old.d.clinical.sources.old.http(base+'/responses',request,env['OPENAI_API_KEY'],timeout=300))
            save(raw,{'request_sha256':sha,'elapsed_seconds':time.monotonic()-start,'response':data})
            return data
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as exc:
            code=getattr(exc,'code',None)
            save(folder/f'transport_error_{tag}_{time.time_ns()}.json',{'status':code,'error_type':type(exc).__name__,
                'attempt':attempt+1,'request_sha256':sha,'billing_may_be_unrecorded':True})
            if attempt==2 or code is not None and code not in [408,429,500,502,503,504,520,522,524]:raise
            time.sleep(5*(attempt+1))


def infer(folder,prompt,packet,schema,validator,env,budget):
    if (folder/'result.json').exists():return
    errors=[]
    for repair in range(2):
        instructions=prompt
        if errors:instructions+='\nCorrect only these format/source-ID errors; do not change the evidentiary task:\n'+json.dumps(errors)
        req={'model':MODEL,'reasoning':{'effort':'medium'},'store':False,'max_output_tokens':budget,
            'instructions':instructions,'input':json.dumps(packet,ensure_ascii=False),
            'text':{'format':{'type':'json_schema','name':'reason_reflection','strict':True,'schema':schema}}}
        tag=f'r{repair}';response=api(folder,tag,req,env);selected=tag
        if response.get('status')=='incomplete' and (response.get('incomplete_details') or {}).get('reason')=='max_output_tokens':
            larger={**req,'max_output_tokens':max(10000,budget*2)};selected=tag+'_capacity'
            response=api(folder,selected,larger,env)
        require(response.get('status')=='completed','Incomplete/refused response is not a negative decision')
        output=json.loads(''.join(c.get('text','') for i in response.get('output',[]) for c in i.get('content',[]) if c.get('type')=='output_text'))
        try:validator(output,packet)
        except (ValueError,KeyError,TypeError) as exc:
            errors=[str(exc)];save(folder/f'validation_{repair}.json',{'valid':False,'errors':errors})
            if repair==1:raise
            continue
        save(folder/f'validation_{repair}.json',{'valid':True})
        save(folder/'result.json',{'output':output,'selected_response_tag':selected,'request_sha256':digest(old.d.clinical.sources.old.encoded(
            load(folder/(selected+'_request.json'))))})
        return


def validate_saved(folder,validator):
    packet=load(folder/'input.json');record=load(folder/'result.json');tag=record['selected_response_tag']
    raw=load(folder/(tag+'_raw_response.json'));request=load(folder/(tag+'_request.json'))
    require(raw['response'].get('status')=='completed','Selected response incomplete')
    require(raw['request_sha256']==record['request_sha256']==digest(old.d.clinical.sources.old.encoded(request)),'Request hash mismatch')
    require(json.loads(request['input'])==packet,'Request used different clinical input')
    require(request['text']['format']['schema']==load(folder/'schema.json'),'Different output schema')
    require(request['model']==MODEL,'Different model')
    actual=json.loads(''.join(c.get('text','') for i in raw['response'].get('output',[]) for c in i.get('content',[]) if c.get('type')=='output_text'))
    require(actual==record['output'],'Stored judgement differs from API response')
    validator(actual,packet)
    return actual


def dispatch(rows,worker,stage):
    env=old.d.clinical.environment();require(bool(env.get('OPENAI_API_KEY')),'Missing credential')
    pending={};errors=[];queue=iter(rows);exhausted=False
    def submit(pool):
        nonlocal exhausted
        if exhausted or errors or (RUN/('STOP_'+stage)).exists():return False
        try:row=next(queue)
        except StopIteration:exhausted=True;return False
        pending[pool.submit(worker,row,env)]=row;return True
    with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
        for _ in range(WORKERS):submit(pool)
        while pending:
            done,_=concurrent.futures.wait(pending,return_when=concurrent.futures.FIRST_COMPLETED)
            for task in done:
                row=pending.pop(task)
                try:task.result()
                except Exception as exc:
                    errors.append({'type':type(exc).__name__,'message':str(exc)[:250],'token':row.get('token')})
                    print(json.dumps({'technical_error':errors[-1]}),flush=True)
            while len(pending)<WORKERS and submit(pool):pass
    require(not errors,'All active calls saved; stopped on technical errors '+json.dumps(errors))
    require(exhausted,'Paused new submissions; active responses retained')
