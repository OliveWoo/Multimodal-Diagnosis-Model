"""Prevent network access during offline validation (not production runtime)."""
import os
import socket
os.environ['PYTHONDONTWRITEBYTECODE']='1'
for key in list(os.environ):
    if key in {'OPENAI_API_KEY','NCBI_API_KEY','PUBMED_API_KEY','ANTHROPIC_API_KEY'}:
        os.environ.pop(key,None)

def blocked(*args,**kwargs):
    raise RuntimeError('Network access is disabled by the repository offline test guard')

socket.socket.connect=blocked
socket.socket.connect_ex=blocked
socket.create_connection=blocked

