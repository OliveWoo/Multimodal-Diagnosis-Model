"""Recover interrupted HTTP bodies, without changing frozen model requests.

Start only after the original dispatcher drains. Does not patch frozen sources.
"""
import http.client
import json
import ssl
import time
from pathlib import Path
import reason_common as c
import reason_judge as j

original_api=c.api


def body_resilient_api(folder,tag,request,env):
    for attempt in range(3):
        try:return original_api(folder,tag,request,env)
        except (http.client.HTTPException,ConnectionError,ssl.SSLError,json.JSONDecodeError) as exc:
            stamp=str(time.time_ns());partial=getattr(exc,'partial',None)
            if isinstance(partial,bytes):c.write_new(folder/('incomplete_http_body_'+stamp+'.bin'),partial)
            c.save(folder/('http_body_recovery_'+stamp+'.json'),{
                'error_type':type(exc).__name__,'attempt':attempt+1,
                'request_sha256':c.digest(c.old.d.clinical.sources.old.encoded(request)),
                'clinical_request_unchanged':True,'billing_for_interrupted_response_unknown':True})
            if attempt==2:raise
            time.sleep(5*(attempt+1))


if __name__=='__main__':
    c.verify('judge_inputs')
    c.save(c.RUN/'judge_transport_amendment_01.json',{
        'cause':'One initial IncompleteRead after 121 completed judgements; original dispatcher drained before restart.',
        'recover_only_interrupted_http_bodies_or_connections':True,
        'same_frozen_prompt_inputs_schema_model_and_selector':True,
        'additional_attempts_max':3,'backoff_seconds':[5,10],
        'script_sha256':c.digest(Path(__file__).read_bytes()),
        'unrecorded_billing_possible_for_initial_interrupted_response':True})
    c.api=body_resilient_api
    j.run()
