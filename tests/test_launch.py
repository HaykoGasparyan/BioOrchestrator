from types import SimpleNamespace
from pathlib import Path
import pytest
from biocentric.launch import read_key, verify_connection, child_environment, LaunchError
from biocentric.service import Service


def test_key_prompt_fails_closed_without_private_terminal():
    with pytest.raises(LaunchError, match='interactive Terminal'):
        read_key({}, interactive=False, prompt=lambda _: pytest.fail('must not prompt'))
    assert read_key({'OPENAI_API_KEY':'test-secret'}, interactive=False)=='test-secret'
    assert read_key({}, interactive=True, prompt=lambda _: 'test-secret')=='test-secret'


def test_connection_is_bounded_and_errors_do_not_expose_key():
    requests=[]
    class Client:
        def __init__(self, **kwargs): requests.append(kwargs)
        def __enter__(self): return self
        def __exit__(self, *args): pass
        @property
        def responses(self): return self
        def create(self, **kwargs):
            requests.append(kwargs)
            return SimpleNamespace(status='completed', output_text='CONNECTED')
    assert verify_connection('test-secret','gpt-4.1-mini',Client)['status']=='connected'
    assert requests[0]['max_retries']==0 and requests[0]['base_url']=='https://api.openai.com/v1'
    assert requests[1]['max_output_tokens']==16 and requests[1]['store'] is False
    class Failing(Client):
        def create(self, **kwargs): raise RuntimeError('secret-body test-secret')
    with pytest.raises(LaunchError) as failure: verify_connection('test-secret','gpt-4.1-mini',Failing)
    assert 'test-secret' not in str(failure.value) and 'secret-body' not in str(failure.value)


def test_child_environment_is_private_and_official(monkeypatch,tmp_path):
    monkeypatch.delenv('OPENAI_API_KEY',raising=False)
    monkeypatch.setenv('OPENAI_BASE_URL','https://example.invalid')
    environment=child_environment('test-secret','gpt-4.1-mini',tmp_path)
    import os
    assert 'OPENAI_API_KEY' not in os.environ
    assert environment['OPENAI_API_KEY']=='test-secret'
    assert environment['OPENAI_BASE_URL']=='https://api.openai.com/v1'
    assert not list(tmp_path.iterdir())


def test_public_example_enables_model_without_changing_private_default(tmp_path):
    service=Service(tmp_path)
    public=service.create_example(allow_model=True,model='configured-model')
    private=service.create_example()
    spec=service.snapshot(public)['project']['spec']
    assert spec['allow_model'] and not spec['local_only'] and spec['model']=='configured-model'
    assert spec['max_model_calls']==3 and spec['max_model_output_tokens']==1500
    assert service.snapshot(private)['project']['spec']['local_only']
    assert not service.snapshot(private)['project']['spec']['allow_model']
    assert service.snapshot(public)['jobs']==[]
