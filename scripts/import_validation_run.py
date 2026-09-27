"""Import verified local batch outputs as a clearly labeled, read-only app run.

This imports existing calculations; it never invokes a model or docking engine.
The source run remains intact and repeated imports reuse the same project.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import io
import json
import math
import os
from pathlib import Path
import re
import uuid

from biocentric import science
from biocentric.models import ProjectSpec
from biocentric.service import Service
from biocentric.store import digest, dump
from biocentric.worker import Pipeline


def _load(content, label):
    try:
        result = json.loads(content)
    except (ValueError, UnicodeDecodeError) as error:
        raise ValueError(f'{label} is not a valid JSON document') from error
    if not isinstance(result, dict):
        raise ValueError(f'{label} must be a JSON object')
    return result


def _checked_bytes(path, inventory):
    path = Path(path).resolve(strict=True)
    record = inventory.get(str(path))
    if record is None:
        raise ValueError(f'Supplied source is missing from the report evidence inventory: {path.name}')
    content = path.read_bytes()
    if digest(content) != record['sha256'] or len(content) != record['bytes']:
        raise ValueError(f'Source artifact integrity failed: {path.name}')
    return content


def _molecules(content, label):
    from rdkit import Chem
    molecules = list(Chem.ForwardSDMolSupplier(io.BytesIO(content), removeHs=False))
    if not molecules or any(mol is None for mol in molecules):
        raise ValueError(f'{label} contains no complete valid molecular records')
    return molecules


def _validate_sources(report_path, preparation_path, docking_path, library_path):
    """Validate the complete evidence bundle before creating any app state."""
    report_bytes = report_path.read_bytes()
    report = _load(report_bytes, 'Validation report')
    if report.get('status') != 'completed':
        raise ValueError('Only a completed validation report can be imported as a completed app run')
    if report.get('target') != 'Curated 1IEP mouse ABL pocket; bound ligand defines box only':
        raise ValueError('This importer supports the recorded 1IEP mouse ABL validation baseline only')
    resources = report.get('resources') or {}
    if any(resources.get(key) != 0 for key in ('llm_calls', 'network_requests')) or resources.get('gpu_used') is not False:
        raise ValueError('Source report does not establish the supported local CPU/no-model validation scope')
    budget = resources.get('budget_seconds')
    if not isinstance(budget, int) or isinstance(budget, bool) or not 10 <= budget <= 86400:
        raise ValueError('Source report must preserve its bounded execution budget')
    if report.get('adversarial_passed') is not True or (report.get('checkpoint_reuse') or {}).get('status') != 'passed':
        raise ValueError('Source validation report has incomplete preparation/checkpoint checks')
    base = report_path.parent
    expected_docking = (base / 'baseline' / 'docking.json').resolve()
    if docking_path != expected_docking:
        raise ValueError('Use canonical baseline/docking.json; *_result.json may retain a stale self-manifest after checkpoint reuse')
    inventory = {}
    for record in report.get('evidence_inventory', []):
        path = Path(record['path']).resolve(strict=True)
        if not path.is_relative_to(base) or not path.is_file():
            raise ValueError('Validation inventory must contain files inside its source run directory')
        if not re.fullmatch(r'[a-f0-9]{64}', str(record.get('sha256', ''))) or not isinstance(record.get('bytes'), int) or isinstance(record.get('bytes'), bool):
            raise ValueError('Validation inventory requires exact hashes and byte counts')
        if str(path) in inventory:
            raise ValueError('Validation inventory contains a duplicate source path')
        inventory[str(path)] = dict(record)
    if not inventory:
        raise ValueError('Validation report has no evidence inventory')
    for path in inventory:
        _checked_bytes(path, inventory)
    primary_bytes = {name: _checked_bytes(path, inventory) for name, path in
                     [('preparation', preparation_path), ('docking', docking_path), ('library', library_path)]}
    preparation = _load(primary_bytes['preparation'], 'Preparation manifest')
    docking = _load(primary_bytes['docking'], 'Docking manifest')
    if preparation.get('status') != 'completed' or docking.get('status') != 'completed':
        raise ValueError('Preparation and docking must have supported completed status; partial runs cannot be promoted')
    compounds, poses = preparation.get('compounds'), docking.get('poses')
    if not compounds or not poses:
        raise ValueError('Import requires actual prepared compounds and saved docking poses')
    manifests = {}
    for result in (preparation, docking):
        for artifact in result.get('artifacts', []):
            path = str(Path(artifact['path']).resolve(strict=True))
            record = inventory.get(path)
            if not record or record['sha256'] != artifact.get('sha256'):
                raise ValueError('Source manifest artifact differs from the verified report inventory')
            manifests[(path, artifact['sha256'])] = record
    library_sha = digest(primary_bytes['library'])
    if any(row.get('source_sha256') != library_sha for row in compounds):
        raise ValueError('Selected input library differs from preparation provenance')
    prepared_path = Path(preparation['prepared_sdf']).resolve(strict=True)
    prepared_bytes = _checked_bytes(prepared_path, inventory)
    prepared_sha = digest(prepared_bytes)
    if (str(prepared_path), prepared_sha) not in manifests:
        raise ValueError('Prepared library is missing from its source manifest')
    identities = {}
    for compound in compounds:
        canonical = compound.get('canonical_smiles')
        if not isinstance(canonical, str) or compound.get('compound_id') != 'cmp_' + digest(canonical.encode())[:16]:
            raise ValueError('Prepared compound identifier does not match its recorded chemical identity')
        if compound['compound_id'] in identities:
            raise ValueError('Prepared compound identities must be unique')
        identities[compound['compound_id']] = canonical
    if Counter(science._canonical(mol) for mol in _molecules(prepared_bytes, 'Prepared library')) != Counter(identities.values()):
        raise ValueError('Prepared SDF identities differ from the preparation manifest')

    request_path = (base / 'baseline' / 'request.json').resolve(strict=True)
    request = _load(_checked_bytes(request_path, inventory), 'Native request')
    request_identity = {key: value for key, value in request.items() if key not in {'fingerprint', 'box'}}
    fingerprint = digest(json.dumps(request_identity, sort_keys=True).encode())
    metadata, baseline = docking.get('metadata') or {}, (report.get('protocols') or {}).get('baseline') or {}
    if not re.fullmatch(r'[a-f0-9]{64}', fingerprint) or any(value != fingerprint for value in
            (request.get('fingerprint'), metadata.get('fingerprint'), baseline.get('protocol_fingerprint'))):
        raise ValueError('Native request, docking and source report do not bind one protocol fingerprint')
    if baseline.get('status') != 'completed' or request.get('library') != prepared_sha:
        raise ValueError('Source baseline is incomplete or uses a different prepared library')
    config = metadata.get('config') or {}
    if request.get('config') != config or metadata.get('box') != request.get('box'):
        raise ValueError('Docking metadata differs from the recorded request configuration')
    declared_config = baseline.get('config') or {}
    if not declared_config or any(config.get(key) != value for key, value in declared_config.items()):
        raise ValueError('Source report configuration differs from the recorded baseline')
    if config.get('cpu') != resources.get('cpu_per_vina_process') or config.get('require_posebusters') is not True:
        raise ValueError('Recorded CPU or physical-check requirements differ from validation scope')
    pinned = report.get('pinned_inputs') or {}
    for role, key in [('receptor', 'receptor_sha256'), ('reference', 'reference_sha256'),
                      ('receptor_pdb', 'receptor_pdb_sha256'), ('adapter_sha256', 'science_adapter_sha256')]:
        if not request.get(role) or request.get(role) != pinned.get(key):
            raise ValueError('Native input hashes differ from the source report pinned inputs')
    if metadata.get('vina') != request.get('vina') or metadata.get('versions') != request.get('versions'):
        raise ValueError('Tool version metadata differs from the source native request')

    completed_ids, pose_ids = set(), set()
    scores = []
    for pose in poses:
        compound_id = pose.get('compound_id')
        if compound_id not in identities or pose.get('canonical_smiles') != identities[compound_id]:
            raise ValueError('Docking pose contains an unrequested or changed compound identity')
        pose_id = pose.get('pose_id')
        if not pose_id or pose_id in pose_ids:
            raise ValueError('Docking poses need unique recorded pose identifiers')
        pose_ids.add(pose_id)
        if (pose.get('engine') not in {'AutoDock Vina', 'vina'} or pose.get('engine_version') != request.get('vina')
                or pose.get('score_units') != 'kcal/mol' or pose.get('score_direction') != 'lower_is_better'):
            raise ValueError('Only consistent native Vina score semantics are supported by this importer')
        score = pose.get('score')
        if not isinstance(score, (int, float)) or isinstance(score, bool) or not math.isfinite(score):
            raise ValueError('Docking scores must be finite native numeric values')
        if pose.get('protocol_fingerprint') != fingerprint or pose.get('receptor_sha256') != request['receptor']:
            raise ValueError('A pose belongs to a different receptor or protocol')
        path = str(Path(pose['pose_path']).resolve(strict=True))
        if (path, pose.get('pose_sha256')) not in manifests:
            raise ValueError('A docking pose is missing from the verified source manifest')
        molecules = _molecules(_checked_bytes(path, inventory), 'Docking pose')
        if len(molecules) != 1 or science._canonical(molecules[0]) != identities[compound_id]:
            raise ValueError('Docking pose SDF identity differs from its prepared compound')
        molecule = molecules[0]
        if not molecule.HasProp('BIOCENTRIC_VINA_SCORE_KCAL_MOL') or not math.isclose(molecule.GetDoubleProp('BIOCENTRIC_VINA_SCORE_KCAL_MOL'), score, abs_tol=1e-8, rel_tol=1e-10):
            raise ValueError('Docking score differs from its registered native pose SDF')
        validity = pose.get('validation') or {}
        validity_artifact = validity.get('artifact') or {}
        validation_path = str(Path(validity_artifact.get('path', '')).resolve())
        if (validation_path, validity_artifact.get('sha256')) not in manifests:
            raise ValueError('Physical validity report is missing from verified source evidence')
        physical = _load(_checked_bytes(validation_path, inventory), 'Physical validity report')
        if {key: value for key, value in validity.items() if key != 'artifact'} != physical:
            raise ValueError('Pose physical-validity flags differ from the source check report')
        checks = physical.get('checks') or {}
        if (physical.get('status') != 'completed' or physical.get('tool') != 'PoseBusters'
                or not checks or any(not isinstance(value, bool) for value in checks.values())
                or physical.get('passed') is not all(checks.values())):
            raise ValueError('Physical validity check report is incomplete or internally inconsistent')
        completed_ids.add(compound_id)
        scores.append(score)
    failures = docking.get('failures') or []
    failed_ids = [row.get('compound_id') for row in failures]
    if len(set(failed_ids)) != len(failed_ids) or any(not row.get('reason') or not row.get('code') for row in failures):
        raise ValueError('Calculation failures need unique compound identifiers and explicit reasons')
    if completed_ids & set(failed_ids) or completed_ids | set(failed_ids) != set(identities):
        raise ValueError('Every requested compound must have completed poses or an explicit failure; no partial import')
    valid = sum(pose['validation']['passed'] is True for pose in poses)
    expected = {'requested_ligands': len(identities), 'ligands_with_poses': len(completed_ids), 'pose_count': len(poses),
                'valid_poses': valid, 'invalid_or_unverified_poses': len(poses) - valid,
                'failed_ligands': len(failed_ids), 'not_completed_ligands': 0,
                'score_units': 'kcal/mol', 'score_direction': 'lower_is_better', 'score_range': [min(scores), max(scores)]}
    if any(baseline.get(key) != value for key, value in expected.items()):
        raise ValueError('Source report counts or score range differ from the verified baseline results')
    preparation_summary = report.get('preparation') or {}
    if (preparation_summary.get('accepted') != len(compounds)
            or preparation_summary.get('excluded') != len(preparation.get('exclusions', []))
            or preparation_summary.get('requested') != preparation.get('input_records')
            or pinned.get('count') != preparation.get('input_records')
            or (report.get('selection') or {}).get('requested') != preparation.get('input_records')):
        raise ValueError('Preparation counts differ from the source report')
    source_binding = {'schema_version': 1, 'source_report_sha256': digest(report_bytes),
                      'source_preparation_sha256': digest(primary_bytes['preparation']),
                      'source_docking_sha256': digest(primary_bytes['docking']), 'source_library_sha256': library_sha,
                      'protocol_fingerprint': fingerprint}
    return {'report': report, 'report_bytes': report_bytes, 'preparation': preparation, 'docking': docking,
            'inventory': inventory, 'prepared_path': str(prepared_path), 'prepared_sha256': prepared_sha,
            'request': request, 'source_binding': source_binding, 'bundle_sha256': digest(dump(source_binding).encode())}


def import_run(report_path, preparation_path, docking_path, library_path, root=None):
    report_path, preparation_path, docking_path, library_path = [Path(path).resolve(strict=True) for path in
                                                               (report_path, preparation_path, docking_path, library_path)]
    source = _validate_sources(report_path, preparation_path, docking_path, library_path)
    report, preparation, docking = source['report'], source['preparation'], source['docking']
    report_hash, bundle_hash = source['source_binding']['source_report_sha256'], source['bundle_sha256']

    service = Service(root)
    for project in service.list_projects():
        if project.get('data', {}).get('research', {}).get('validation_run', {}).get('source_report_sha256') == report_hash:
            previous = project['data']['research']['validation_run']
            if previous.get('source_bundle_sha256') != bundle_hash:
                raise ValueError('This report already has an import with different source bindings; inspect it before retrying')
            if project['state'] == 'completed':
                for artifact in service.store.artifacts(project['id']):
                    service.store.read_artifact(project['id'], artifact)
                return project['id']
            raise ValueError('An incomplete import of this report already exists; inspect it before retrying')
    count = len(preparation['compounds'])
    spec = ProjectSpec(name=f'Local validation · {count} diverse ligands', target='P00520', organism='Mus musculus',
                       context='Recorded software validation on chemically diverse compounds. No LLM, literature-driven tool choice, measured binding, or comparative engine benchmark is claimed.',
                       local_only=True, allow_model=False, max_model_calls=0, workflow_mode='research',
                       shortlist_size=10, cpu=source['request']['config']['cpu'], seed=source['request']['config']['seed'],
                       compute_budget_seconds=report['resources']['budget_seconds'])
    steps = [
        {'id': '01', 'label': 'Select the diverse software-test library'},
        {'id': '04', 'label': 'Prepare compounds and preserve exclusions'},
        {'id': '08', 'label': 'Run the recorded Vina batch and pose checks'},
        {'id': '09', 'label': 'Rank results and verify structure downloads'},
    ]
    pid = service.store.create(spec.model_dump(), steps)
    jid = service.store.enqueue(pid)
    owner = uuid.uuid4().hex
    if not service.store.claim(jid, owner_token=owner, worker_pid=os.getpid(), project_id=pid):
        raise ValueError('Import could not acquire its local job')
    service.store.research_update(pid, validation_run={**source['source_binding'], 'source_bundle_sha256': bundle_hash,
        'source_report': str(report_path), 'scope': 'software_validation', 'model_calls': 0,
        'source_preparation': str(preparation_path), 'source_docking': str(docking_path),
        'prepared_compounds': count, 'imported_calculations': True, 'read_only': True,
        'source_resources': report['resources'], 'importer_sha256': science.sha256(Path(__file__))},
        goal=spec.context, questions=[], phase='importing', hardware=report['resources'].get('hardware', 'Not reported'))
    try:
        saved_report = service.store.artifact(pid, 'diverse-validation-report.json', source['report_bytes'], 'validation_report', jid)
        registered = {}
        for path, record in source['inventory'].items():
            kind = 'library' if path == str(library_path) else 'scientific_output'
            registered[(path, record['sha256'])] = service.store.artifact(pid, Path(path).name, _checked_bytes(path, source['inventory']), kind, jid)
        library = registered[(str(library_path), source['source_binding']['source_library_sha256'])]
        prepared_artifact = registered[(source['prepared_path'], source['prepared_sha256'])]
        source_map = [{'source_path': path, 'sha256': sha, 'artifact_id': artifact['id']} for (path, sha), artifact in registered.items()]
        binding_artifact = service.store.artifact(pid, 'import-source-bindings.json', dump({'binding': source['source_binding'], 'files': source_map}).encode(), 'validation_report', jid)
        evidence_ids = {name: next((artifact['id'] for (path, sha), artifact in registered.items()
                                    if Path(path) == report_path.parent / name), None)
                        for name in ('events.jsonl', 'selection.json', 'parameter_comparisons.json',
                                     'checkpoint_reuse.json', 'adversarial_verification.json')}
        prepared = deepcopy(preparation)
        prepared['prepared_artifact_id'] = prepared_artifact['id']
        prepared['prepared_sdf'] = prepared_artifact['path']
        for pose in docking['poses']:
            artifact = registered[(str(Path(pose['pose_path']).resolve()), pose['pose_sha256'])]
            pose['source_provenance'] = {'pose_path': pose['pose_path'], 'job_id': pose.get('job_id'),
                                         'artifact_ids': pose.get('artifact_ids', []), 'pose_artifact_id': pose.get('pose_artifact_id'),
                                         'docking_manifest_sha256': source['source_binding']['source_docking_sha256']}
            pose.update(artifact_ids=[artifact['id']], pose_artifact_id=artifact['id'], pose_path=artifact['path'], job_id=jid)
        service.store.update(pid, data={'preparation': prepared, 'screening': docking,
            'result_status': 'exploratory_software_validation', 'limitations': [
                'NCI diversity samples have no matched target-activity labels in this run.',
                'Only Vina was executed; changing its settings is not a comparison of different docking engines.',
                'A passing software run and pose checks do not establish measured binding.'],
            'source_validation': {'binding_artifact_id': binding_artifact['id'], 'source_report_artifact_id': saved_report['id'],
                                  'target_statement': report['target'], 'protocol': source['request'],
                                  'evidence_artifact_ids': evidence_ids}})

        def completed(sid, reason, output):
            artifact = service.store.artifact(pid, f'validation-stage-{sid}.json', dump(output).encode(), 'research_stage', jid)
            service.store.step(pid, sid, 'completed', reason=reason, outputs={**output, 'result_artifact': artifact['id']})
            service.store.event(pid, kind='stage', stage=sid, name=next(s['label'] for s in steps if s['id'] == sid),
                                status='completed', message=reason, outputs={'result_artifact': artifact['id']},
                                provenance='Imported verified local calculation; not a model response')

        completed('01', f'{count} prepared compounds from a reproducible diverse library; this is a software test, not an assay benchmark.',
                  {'summary': f'Diverse input library: {count} prepared compounds.', 'source_report_artifact_id': saved_report['id'],
                   'library_artifact_id': library['id'], 'validation_report': report})
        completed('04', f'RDKit prepared {count} compounds and retained {len(prepared.get("exclusions", []))} exclusions.',
                  {'summary': 'Actual chemical preparation, including recorded stereochemistry and charge assumptions.',
                   'prepared_count': count, 'exclusions': prepared.get('exclusions', []), 'prepared_artifact_id': prepared_artifact['id']})
        valid = sum(pose.get('validation', {}).get('passed') is True for pose in docking['poses'])
        completed('08', f'Actual Vina results: {len(docking["poses"])} poses, {valid} passing physical checks, {len(docking.get("failures", []))} calculation failures.',
                  {'summary': 'Fixed local Vina software-validation protocol; no AI method choice or predictive benchmark.',
                   'pose_count': len(docking['poses']), 'valid_pose_count': valid, 'failures': docking.get('failures', []),
                   'metadata': docking.get('metadata', {}), 'event_log_artifact_id': evidence_ids['events.jsonl']})
        service.store.event(pid, kind='tool_call', stage='08', name='Imported AutoDock Vina batch', status='completed',
                            inputs=docking.get('metadata', {}).get('config'),
                            outputs={'poses': len(docking['poses']), 'physical_checks_passed': valid, 'source_report': saved_report['id']})
        pipeline = Pipeline(service.store, pid, jid, owner_token=owner)
        ranking = pipeline.prioritization()
        outputs = pipeline.handoff()
        completed('09', f'Ranked {ranking["count"]} ligands; retained {ranking["selected"]} diverse, physically valid poses for software-test review.',
                  {'summary': 'Scores and selected structures are downloadable. These selections are not confirmed binders.',
                   'ranked_count': ranking['count'], 'kept_count': ranking['selected'],
                   'shortlist': outputs['shortlist'], 'limitations': service.store.project(pid)['data']['limitations']})
        service.store.research_update(pid, phase='completed')
        service.store.finish(jid, 'completed', owner_token=owner, project_id=pid)
        service.export(pid)
        for artifact in service.store.artifacts(pid):
            service.store.read_artifact(pid, artifact)
        return pid
    except Exception as error:
        service.store.research_update(pid, phase='import_failed')
        service.store.finish(jid, 'failed', str(error), project_id=pid)
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('report', 'preparation', 'docking', 'library'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--root', default=None)
    args = parser.parse_args()
    print(import_run(args.report, args.preparation, args.docking, args.library, args.root))
