"""Readable summaries of recorded workflow state; these never call a model."""
from __future__ import annotations


def model_activity(snapshot: dict, *, key_available: bool) -> dict:
    spec = snapshot.get('project', {}).get('spec', {})
    data = snapshot.get('data', {})
    steps = snapshot.get('steps', [])
    strategy = next((step for step in steps if step.get('id') == '03'), {})
    advice = (data.get('model_advice') or (strategy.get('outputs') or {}).get('model_advice')
              or (data.get('plan') or {}).get('model_advice') or {})
    events = [event for event in snapshot.get('events', []) if event.get('kind') == 'model_call']
    requests = [event for event in events if event.get('status') == 'running']
    responses = {event['outputs']['response_id']: event['outputs'] for event in events
                 if isinstance(event.get('outputs'), dict) and event['outputs'].get('response_id')}
    usage = [response['usage'] for response in responses.values() if isinstance(response.get('usage'), dict)]
    calls = len(requests)
    rejected_advice = {}
    state, reason = 'waiting', 'The advisor runs at Stage 03, after target evidence is collected.'
    if spec.get('mode') == 'offline_fixture':
        state, reason = 'disabled', 'This is an offline demonstration. It makes no model calls.'
    elif spec.get('local_only', True):
        state, reason = 'disabled', 'This project keeps inputs local, so remote LLM calls are disabled.'
    elif not spec.get('allow_model', False):
        state, reason = 'disabled', 'Language-model use is disabled in this project’s saved settings.'
    elif spec.get('max_model_calls', 3) == 0:
        state, reason = 'disabled', 'This project has a model-call budget of zero.'
    elif advice.get('summary'):
        state, reason = 'completed', 'The advisor returned a structured proposal. Its explanation appears below.'
    elif (data.get('model_repair') or {}).get('status') == 'running':
        state, reason = 'running', 'Correcting the rejected advice with one remaining authorized call.'
    elif advice.get('status') == 'failed' or (events and events[-1].get('status') == 'failed'):
        state, reason = 'failed', advice.get('reason') or events[-1].get('message') or 'The model attempt failed.'
    elif strategy.get('status') == 'running' and calls:
        state, reason = 'running', 'Waiting for model output or reading evidence requested by the model.'
    elif data.get('model_attempted'):
        state, reason = 'incomplete', 'A previous model attempt has no validated final advice. Resume does not make an automatic paid retry.'
    elif not key_available:
        state, reason = 'unconfigured', 'This app process has no API key. Launch it with the private key prompt before running this project.'
    if state == 'failed':
        from .agent import Advice
        for event in reversed(events):
            found = False
            for item in (event.get('outputs') or {}).get('output', []):
                for part in item.get('content', []) if item.get('type') == 'message' else []:
                    if part.get('type') != 'output_text':
                        continue
                    found = True
                    try:
                        rejected_advice = Advice.model_validate_json(part['text']).model_dump()
                    except ValueError:
                        pass  # Malformed public text remains inspectable in its original event.
            if found:
                break
    return {'state': state, 'reason': reason, 'model': spec.get('model', 'gpt-4.1-mini'),
            'calls': calls, 'budget': spec.get('max_model_calls', 3), 'advice': advice, 'rejected_advice': rejected_advice,
            'input_tokens': sum(item.get('input_tokens', 0) or 0 for item in usage) if usage else None,
            'output_tokens': sum(item.get('output_tokens', 0) or 0 for item in usage) if usage else None,
            'key_available': key_available}


def stage_explanation(snapshot: dict, step: dict) -> dict:
    """Explain recorded results and declared rules, without inventing model thoughts."""
    sid, status = step.get('id'), step.get('status', 'pending')
    output = step.get('outputs') or {}
    data = snapshot.get('data', {})
    reason = step.get('reason') or 'No recorded outcome yet.'
    rule = step.get('postprocessing') or 'See the recorded stage contract.'
    details = []
    if status in {'completed', 'complete', 'succeeded'}:
        if sid == '01' and output.get('identity'):
            identity = output['identity']
            reason = f"Used target {identity.get('accession') or identity.get('name', 'provided sequence')} in {identity.get('organism', 'the recorded organism')}."
            details.append(identity.get('verification') or 'Identity source and mapping are recorded in the output.')
        elif sid == '02':
            reason = f"Retained {output.get('count', 0)} assay records with their original measurement context."
            details.append(output.get('control_status', 'Control availability is recorded in the output.'))
        elif sid == '03':
            reason = output.get('rationale') or (data.get('plan') or {}).get('rationale') or reason
            details.append(f"Workflow route: {output.get('strategy', 'not recorded').replace('_', ' ')}.")
        elif sid == '04':
            reason = 'The supplied structures passed the supported-component checks and the recorded protocol review.'
            if (output.get('review') or {}).get('ph_assumption') is not None:
                details.append(f"Recorded pH assumption: {output['review']['ph_assumption']}.")
        elif sid == '05':
            reason = f"Prepared {len(output.get('compounds', []))} compound records; {len(output.get('exclusions', []))} exclusions were recorded."
        elif sid == '06':
            reason = f"The pilot produced {len(output.get('poses', []))} poses and recorded {len(output.get('failures', []))} failures."
        elif sid == '07':
            reason = 'Screening remains exploratory because no independent held-out benchmark established predictive performance.'
        elif sid == '08':
            reason = f"Screening produced {output.get('count', 0)} poses and recorded {len(output.get('failures', []))} failures."
        elif sid == '09':
            reason = f"Selected {output.get('selected', 0)} of {output.get('count', 0)} candidates using pose validity, Vina rank, and fingerprint diversity."
        elif sid == '10':
            reason = 'Exported the computational shortlist, provenance, and proposed experimental follow-up.'
            details.append(output.get('proposal', 'Experimental work still requires scientist review.'))
        details.extend(output.get('limitations') or [])
    if sid == '03':
        actor = 'Workflow rules select the route; GPT can add an advisory proposal when enabled.'
    else:
        actor = 'Scientific tools and workflow rules perform this stage. No LLM is called here.'
    return {'actor': actor, 'reason': reason, 'rule': rule, 'details': details}
