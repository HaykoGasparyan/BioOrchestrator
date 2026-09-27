"""Reviewed, bounded literature ranking benchmarks for the discovery workspace.

Experimental p-values remain endpoint/assay specific. Native scores are never
relabelled as affinity. A separate process owns each immutable accepted run.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import uuid

from . import discovery_evidence as evidence, science
from . import pocket_identification as pockets

FACTORS = {'M': 1., 'mM': 1e-3, 'uM': 1e-6, 'µM': 1e-6, 'μM': 1e-6, 'nM': 1e-9, 'pM': 1e-12}
SCORES = {'vina_score': ('Vina score', 'kcal/mol', False),
          'gnina_cnn': ('CNNscore', 'dimensionless', True),
          'gnina_affinity': ('CNNaffinity', 'pK (model)', True),
          'biocentric_affinity': ('BioCentric affinity', 'pK (model)', True)}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def runtime_fingerprint(candidates):
    from .engines import _executable
    model_files = {}
    tools = {p['docking']['tool'] for p in candidates}
    tools.update('vina' if p['scoring']['tool'] == 'vina_score' else 'gnina'
                 for p in candidates if p['scoring']['tool'] != 'biocentric_affinity')
    paths = [Path(__file__), Path(science.__file__), Path(__file__).with_name('engines.py'), Path(__file__).with_name('engine_workflow.py')]
    paths.append(Path(pockets.__file__))
    paths.extend(Path(__file__).with_name(name) for name in ('bayesian_search.py', 'discovery_search_runner.py', 'pipeline_decisions.py', 'evaluation_archive.py'))
    if any(p['scoring']['tool'] == 'biocentric_affinity' for p in candidates):
        from . import biocentric_model, discovery_affinity
        runtime = biocentric_model.installation()
        if runtime['blockers']:
            raise ValueError('; '.join(runtime['blockers']))
        model_files = {str(Path(item['path']).resolve()): item['sha256'] for item in runtime['files'].values()}
        paths.extend([Path(runtime['python']), Path(biocentric_model.__file__),
                      Path(discovery_affinity.__file__), Path(__file__).with_name('biocentric_predictor.py')])
    if any(p.get('pocket', {}).get('tool') == 'fpocket' for p in candidates):
        if not pockets.executable():
            raise ValueError('fpocket runtime unavailable')
        paths.append(Path(pockets.executable()))
    for name in tools:
        executable = _executable(name)
        if not executable:
            raise ValueError(f'{name} runtime is unavailable')
        paths.append(Path(executable))
        manifest = Path(__file__).resolve().parents[1] / '.tools/linux/gpu-runtime-manifest.json'
        if name == 'gnina' and Path(executable) == manifest.parent / 'gnina' and manifest.is_file():
            paths.append(manifest)
            for item in json.loads(manifest.read_text())['files']:
                path = Path(item['path']).resolve()
                if not path.is_relative_to(manifest.parent) or science.runtime_sha256(path) != item['sha256']:
                    raise ValueError('GNINA runtime differs from the verified installation')
                paths.append(path)
    return {'files': {**model_files, **{str(p.resolve()): science.runtime_sha256(p) for p in paths}},
            'versions': {name: science._version(name) for name in ('rdkit', 'meeko', 'posebusters', 'scipy')}}


def normalize(rows, organism='Homo sapiens'):
    """Retain originals and explicit exclusions; do not turn bounds into labels."""
    from rdkit import Chem
    records, excluded = [], []
    seen = set()
    for row in rows:
        raw = row.get('raw') or {}
        source_id = str(raw.get('activity_id') or row.get('source_url') or digest(row))
        if source_id in seen:
            continue
        seen.add(source_id)
        record = {**deepcopy(row), 'source_id': source_id, 'p_value': None, 'p_relation': None, 'eligible': False}
        reasons = []
        try:
            if row.get('endpoint') not in {'Kd', 'Ki', 'IC50'} or row.get('units') not in FACTORS:
                raise ValueError('Unsupported endpoint or concentration units')
            if isinstance(row.get('value'), bool):
                raise ValueError('Invalid measurement')
            molar = float(row['value']) * FACTORS[row['units']]
            if not math.isfinite(molar) or molar <= 0:
                raise ValueError('Measurement must be finite and positive')
            relation = row.get('relation')
            if relation not in {'=', '<', '<=', '>', '>=', '~'}:
                raise ValueError('Missing or unsupported measurement relation')
            record.update(p_value=-math.log10(molar), p_relation={'<': '>', '<=': '>=', '>': '<', '>=': '<='}.get(relation, relation), p_endpoint='p' + row['endpoint'])
            if relation != '=':
                reasons.append('Censored/approximate measurement; excluded from exact ranking')
        except (ValueError, TypeError, OverflowError) as exc:
            reasons.append(str(exc))
        if (row.get('organism') or '').casefold() != organism.casefold():
            reasons.append('Assay species missing or mismatched')
        if row.get('_assay_checked') and (row.get('assay_metadata') or {}).get('confidence_score') != 9:
            reasons.append('Assay-to-target assignment is not confirmed as direct single protein')
        if row.get('validity') not in {None, '', 'Manually validated'}:
            reasons.append('Source quality flag: ' + row['validity'])
        if raw.get('potential_duplicate') in {1, '1', True}:
            reasons.append('Potential duplicate flagged by source')
        if not row.get('assay') or not row.get('compound_id') or not row.get('document'):
            reasons.append('Missing compound, assay or publication provenance')
        mol = Chem.MolFromSmiles(row.get('smiles') or '')
        if mol is None or mol.GetNumAtoms() == 0:
            reasons.append('No valid molecular structure')
        else:
            try:
                science._supported(mol)
                record['canonical_smiles'] = science._canonical(mol)
            except science.ScienceError as exc:
                reasons.append(str(exc))
        # One assay, endpoint, target and variant; missing construct is disclosed.
        record['cohort'] = digest([row.get('target_id'), row.get('assay'), row.get('endpoint'), row.get('organism'), raw.get('assay_variant_accession'), raw.get('assay_variant_mutation')])
        record['reasons'] = reasons
        record['eligible'] = not reasons
        records.append(record)
    groups = defaultdict(list)
    for record in records:
        if record['eligible']:
            groups[(record['cohort'], record['canonical_smiles'])].append(record)
    for duplicates in groups.values():
        if len({r['p_value'] for r in duplicates}) > 1:
            for record in duplicates:
                record['eligible'] = False
                record['reasons'].append('Conflicting replicate values in the same assay')
        else:
            for record in duplicates[1:]:
                record['eligible'] = False
                record['reasons'].append('Duplicate molecular identity in the same assay')
    return records


def fetch_dataset(identity, root, first_page=None, max_records=500, progress=None):
    if isinstance(max_records, bool) or not isinstance(max_records, int) or max_records < 1:
        raise ValueError('The database search maximum must be a positive whole number')
    first = first_page if first_page and first_page.get('page', 0) == 0 and not first_page.get('error') else evidence.binding_page(identity, root, 0)
    rows = list(first.get('rows', []))
    total = first.get('total') or 0
    errors = []
    for page in range(1, math.ceil(min(total, max_records) / 50)):
        if progress:
            progress(min(page * 50, max_records), min(total, max_records))
        try:
            rows.extend(evidence.binding_page(identity, root, page)['rows'])
        except (evidence.RetrievalError, OSError, ValueError) as exc:
            errors.append(f'Page {page + 1}: {exc}')
    rows = deepcopy(rows[:max_records])
    assay_ids = sorted({r['assay'] for r in rows if re.fullmatch(r'CHEMBL[0-9]+', r.get('assay') or '')})
    assays = {}
    for start in range(0, len(assay_ids), 50):
        try:
            result = evidence.client(root)._request('ChEMBL', 'https://www.ebi.ac.uk/chembl/api/data/assay.json',
                {'assay_chembl_id__in': ','.join(assay_ids[start:start + 50]), 'limit': 50})
            assays.update({a['assay_chembl_id']: a for a in result['data'].get('assays', [])})
        except evidence.RetrievalError as exc:
            errors.append('Assay metadata: ' + str(exc))
    for row in rows:
        assay = assays.get(row.get('assay'), {})
        raw = row.get('raw') or {}
        row.update(_assay_checked=True, assay_metadata=assay)
        row['organism'] = assay.get('assay_organism')
        if (not row['organism'] and assay.get('confidence_score') == 9 and assay.get('target_chembl_id') == row.get('target_id')
                and str(raw.get('target_tax_id')) == str(identity.get('taxon_id'))):
            row['organism'] = raw.get('target_organism')
            row['organism_basis'] = 'Direct single-protein target species; assay organism unreported'
    rows = rows[:max_records]
    records = normalize(rows, identity['organism'])
    return {'identity': deepcopy(identity), 'records': records, 'source_total': total,
            'loaded': len(records), 'cap': max_records, 'errors': errors,
            'complete': len(rows) >= total and not errors, 'created_at': science._now(),
            'endpoints': dict(Counter(r.get('endpoint') or 'Unknown' for r in records)),
            'eligible': sum(r['eligible'] for r in records)}


def cohorts(dataset):
    groups = defaultdict(list)
    for r in dataset['records']:
        if r['eligible']:
            groups[r['cohort']].append(r)
    return dict(sorted(groups.items(), key=lambda pair: (-len(pair[1]), pair[0])))


def mutation_set(text):
    text = (text or '').upper().strip()
    if text in {'', 'NO', 'NONE', 'WILD TYPE', 'WILD-TYPE'}:
        return set()
    values = set(re.findall(r'\b[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]\b', text))
    return values if values else None


def assay_mutations(row):
    raw = row.get('raw') or {}
    variant = (row.get('assay_metadata') or {}).get('variant_sequence') or {}
    declared = raw.get('assay_variant_mutation') or variant.get('mutation')
    if declared:
        return mutation_set(declared)
    description = row.get('assay_description') or ''
    found = mutation_set(description)
    return found if found else None if 'mutant' in description.casefold() else set()


def compatible_reference(row, reference_record):
    title = (reference_record.get('title') or '').upper()
    description = (row.get('assay_description') or '').upper()
    domain = next((d for d in ('SH2', 'SH3') if re.search(r'\b' + d + r'\b', title)), None)
    if domain and not re.search(r'\b' + domain + r'\b', description):
        return False
    if 'KINASE DOMAIN' in title and re.search(r'\bSH[23]\b', description):
        return False
    assay = assay_mutations(row)
    structure = mutation_set(', '.join(reference_record.get('mutations', [])))
    return assay is not None and structure is not None and assay == structure


def reference_bundle(record, component, chain, root):
    """Extract deposited ligand coordinates, bind CCD chemistry and target chain."""
    import gemmi
    from rdkit import Chem
    from rdkit.Chem import AllChem
    if chain not in record['target_chains'] or component not in {r['id'] for r in record['compounds']}:
        raise ValueError('Reference must belong to the selected experimental target entry')
    geometry = evidence.coordinates(record['id'], root)
    structure = gemmi.make_structure_from_block(gemmi.cif.read_string(geometry['text']).sole_block())
    model = structure[0]
    target = model.find_chain(chain)
    if target is None:
        raise ValueError('Selected target chain has no coordinates')
    # Only one explicitly reviewed protein chain is retained. No atom repair.
    protein = target.clone()
    for i in reversed(range(len(protein))):
        if protein[i].het_flag != 'A':
            del protein[i]
    atoms = [a for r in protein for a in r if not a.is_hydrogen()]
    if not atoms:
        raise ValueError('No standard protein atoms in the selected chain')
    candidates = []
    for ch in model:
        for residue in ch:
            if residue.name == component:
                distance = min(a.pos.dist(b.pos) for a in residue if not a.is_hydrogen() for b in atoms)
                candidates.append((distance, ch.name, residue))
    if not candidates:
        raise ValueError('No deposited coordinates for this ligand')
    distance, ligand_chain, residue = min(candidates, key=lambda r: r[0])
    if distance > 6:
        raise ValueError('Selected ligand is not near the target chain')
    if any(a.altloc not in {'\x00', ' ', ''} for r in protein for a in r) or any(a.altloc not in {'\x00', ' ', ''} for a in residue):
        raise ValueError('Alternate atom locations need an explicit curated reference')
    # Propose only remote incomplete terminal-tail omission, never pocket repair.
    # The exact omissions are shown before approval and frozen in the plan.
    import meeko
    templates = json.loads((Path(meeko.__file__).parent / 'data/residue_chem_templates.json').read_text())['residue_templates']
    incomplete = {}
    for index, res in enumerate(protein):
        # Histidine templates share these heavy atoms; this does not assign its protonation.
        template_name = 'HIE' if res.name == 'HIS' else res.name
        if template_name not in templates:
            raise ValueError('Nonstandard residue requires a curated receptor: ' + res.name)
        expected = {name for name in templates[template_name]['atom_name'] if not name.startswith('H')}
        missing = expected - {a.name for a in res if not a.is_hydrogen()}
        if missing:
            incomplete[index] = sorted(missing)
    terminal = set()
    for indexes in (range(len(protein)), reversed(range(len(protein)))):
        for index in indexes:
            if index not in incomplete:
                break
            terminal.add(index)
    if set(incomplete) - terminal:
        raise ValueError('Internal missing heavy atoms require a curated receptor; no automatic pocket repair')
    omitted = []
    for index in sorted(terminal, reverse=True):
        res = protein[index]
        separation = min(a.pos.dist(b.pos) for a in res if not a.is_hydrogen() for b in residue if not b.is_hydrogen())
        if separation <= 12:
            raise ValueError('Incomplete terminal residue within 12 Å of ligand; curated receptor required')
        omitted.append({'chain': chain, 'residue': str(res.seqid), 'name': res.name, 'missing_atoms': incomplete[index], 'ligand_distance_angstrom': separation})
        del protein[index]
    def pdb(ch):
        s = gemmi.Structure()
        m = gemmi.Model('1')
        m.add_chain(ch)
        s.add_model(m)
        return s.make_pdb_string()
    receptor_text = pdb(protein)
    ligand = gemmi.Chain(ligand_chain)
    ligand.add_residue(residue.clone())
    deposited = Chem.MolFromPDBBlock(pdb(ligand), removeHs=True)
    ccd = evidence.client(root)._request('RCSB CCD', 'https://data.rcsb.org/rest/v1/core/chemcomp/' + component)
    descriptor = ccd['data'].get('rcsb_chem_comp_descriptor') or {}
    template = Chem.MolFromSmiles(descriptor.get('SMILES_stereo') or descriptor.get('SMILES') or '')
    if template is None or deposited is None or not template.GetNumAtoms():
        raise ValueError('Reference chemistry could not be established from CCD')
    ligand_mol = AllChem.AssignBondOrdersFromTemplate(template, deposited)
    if science._canonical(ligand_mol) != science._canonical(template):
        raise ValueError('Deposited ligand chemistry/stereochemistry differs from CCD')
    science._supported(ligand_mol)
    ligand_mol.SetProp('_Name', component)
    reference_text = Chem.MolToMolBlock(ligand_mol) + '\n$$$$\n'
    return {'pdb_id': record['id'], 'component': component, 'chain': chain,
            'ligand_chain': ligand_chain, 'residue': str(residue.seqid), 'distance_angstrom': distance,
            'record': deepcopy(record), 'source_sha256': geometry['sha256'], 'source_url': geometry['url'],
            'ccd_sha256': digest(ccd['data']), 'ccd_url': 'https://data.rcsb.org/rest/v1/core/chemcomp/' + component,
            'receptor_pdb': receptor_text, 'reference_sdf': reference_text,
            'omitted_terminal_residues': omitted,
            'canonical_smiles': science._canonical(ligand_mol)}


def chain_label(pair):
    method = pair.get('pocket', pockets.REFERENCE)['tool']
    name = {'reference_box': 'Crystal pocket', 'p2rank': 'P2Rank', 'deeppocket': 'DeepPocket'}.get(method, method)
    return ' → '.join([name, pair['docking']['tool'], pair['scoring']['tool']])


def measurement_key(row):
    return row.get('measurement_id') or row.get('compound_id')


def combination_progress(review, state):
    """One unit per fully attempted chain, including validation; unrun rows never count."""
    expected = {measurement_key(r) for r in review['records']} | {'crystal:' + review['reference']['component']}
    grouped = defaultdict(list)
    for row in state.get('rows', []):
        grouped[row['combination']].append(row)
    finished = passed = 0
    statuses = {}
    for pair in review['candidates']:
        rows = grouped[pair['id']]
        complete = (len(rows) == len(expected) and {measurement_key(r) for r in rows} == expected
                    and all(r.get('status') in {'completed', 'failed'} for r in rows))
        successful = complete and all(r.get('status') == 'completed' for r in rows)
        finished += int(complete)
        passed += int(successful)
        statuses[pair['id']] = ('Not selected' if rows and all(r.get('status') == 'not_selected' for r in rows) else
                               'Evaluated' if successful else 'Failed' if complete else
                               'Running' if state.get('status') == 'running' and state.get('active_combination') == pair['id'] else
                               'Incomplete' if rows else 'Pending')
    return {'completed': finished, 'successful': passed,
            'total': (review.get('search') or {}).get('max_trials', len(review['candidates'])),
            'pool_total': len(review['candidates']),
            'selected_total': len(review['candidates']) + len(review.get('excluded', [])),
            'unavailable_total': len(review.get('excluded', [])), 'statuses': statuses}


def combination_caption(review):
    progress = combination_progress(review, {})
    return (f"{progress['selected_total']} combinations · {progress['pool_total']} executable"
            + (f" · Bayesian budget {progress['total']}" if review.get('search') else '')
            + (f" · {progress['unavailable_total']} unavailable" if progress['unavailable_total'] else ''))


def evaluation_records(dataset, reference_record):
    """All eligible retrieved records; separate cohorts and explicit reference exclusions."""
    records, excluded = [], []
    for row in dataset['records']:
        if not row['eligible']:
            excluded.append({'source_id': row.get('source_id'), 'reason': '; '.join(row.get('reasons', []))})
        elif not compatible_reference(row, reference_record):
            excluded.append({'source_id': row.get('source_id'), 'reason': 'Crystal/assay domain or mutation annotations incompatible or unresolved'})
        else:
            records.append(row)
    return sorted(records, key=lambda r: (r['cohort'], r['compound_id'], r['source_id'])), excluded


def candidate_protocols(selected):
    from itertools import product
    from .engines import engine_catalog
    catalog = {r['id']: r for r in engine_catalog()}
    runnable, excluded = [], []
    for pocket, dock, score in product(selected.get('pocket_identification', [pockets.REFERENCE]), selected['docking'], selected['scoring']):
        reason = pockets.blocker(pocket)
        if not reason:
            if dock['tool'] not in {'vina', 'gnina', 'unidock'}:
                reason = 'No automatic benchmark input adapter for this docking tool'
            elif catalog[dock['tool']]['blockers']:
                reason = '; '.join(catalog[dock['tool']]['blockers'])
            elif score['tool'] not in SCORES:
                reason = 'Needs a separate validated scoring/trajectory adapter'
            elif score['tool'].startswith('gnina') and catalog['gnina_affinity']['blockers']:
                reason = '; '.join(catalog['gnina_affinity']['blockers'])
            elif score['tool'] == 'vina_score' and catalog['vina']['blockers']:
                reason = '; '.join(catalog['vina']['blockers'])
            elif score['tool'] == 'gnina_cnn' and score['parameters'].get('cnn_scoring') != 'rescore':
                reason = 'Benchmark supports fixed-pose rescoring; refinement needs a separate protocol'
            elif score['tool'] == 'biocentric_affinity':
                from .discovery_affinity import Parameters
                try:
                    par = Parameters.model_validate(score['parameters'])
                    reason = '; '.join(catalog['biocentric_affinity']['blockers'])
                    if not reason and par.device == 'cuda' and not catalog['biocentric_affinity']['hardware']['gpu_ready']:
                        reason = 'BioCentric CUDA device unavailable; select CPU explicitly or configure the GPU'
                except ValueError as exc:
                    reason = str(exc)
        pair = {'id': digest([pocket, dock, score])[:16], 'pocket': deepcopy(pocket), 'docking': deepcopy(dock), 'scoring': deepcopy(score)}
        (excluded if reason else runnable).append({**pair, **({'reason': reason} if reason else {})})
    return runnable, excluded


def ligand_subset(records, max_compounds=None):
    """Freeze a score-independent subset; retain every assay row for chosen IDs."""
    if max_compounds is not None and (type(max_compounds) is not int or max_compounds < 3):
        raise ValueError('The ligand limit must be an integer of at least three')
    ids = {r['compound_id'] for r in records}
    order = lambda key: (digest(['evaluation-ligand-subset-v1', 42, key]), key)
    chosen = set(ids)
    if max_compounds is not None and len(ids) > max_compounds:
        groups = defaultdict(set)
        for row in records:
            groups[row['cohort']].add(row['compound_id'])
        anchor = sorted(groups, key=lambda key: (-len(groups[key]), key))[0]
        if len(groups[anchor]) < 3:
            raise ValueError('A shorter evaluation needs an assay with at least three distinct ligands')
        # Preserve a minimum rankable assay, then sample all other IDs by a
        # stable seeded hash. Never use measured activities or model results.
        chosen = set(sorted(groups[anchor], key=order)[:3])
        chosen.update(sorted(ids - chosen, key=order)[:max_compounds - len(chosen)])
    kept = [r for r in records if r['compound_id'] in chosen]
    excluded = [{'source_id': r.get('source_id'), 'reason': 'Outside the frozen ligand subset'}
                for r in records if r['compound_id'] not in chosen]
    return kept, excluded, {'max_compounds': max_compounds, 'available_compounds': len(ids),
                           'selected_compounds': len(chosen), 'seed': 42,
                           'method': 'Three hash-selected IDs from the largest assay, then seeded hash selection across remaining IDs; all assay records retained for selected IDs. No activity or prediction based selection.'}


def make_review(dataset, cohort_id, selected, reference, max_compounds=None):
    records, record_exclusions = evaluation_records(dataset, reference['record'])
    records, subset_exclusions, subset = ligand_subset(records, max_compounds)
    record_exclusions.extend(subset_exclusions)
    if len(records) < 3 or not any(len(rows) >= 3 for rows in cohorts({'records': records}).values()):
        raise ValueError('At least one compatible assay needs three eligible measurements for ranking')
    candidates, excluded = candidate_protocols(selected)
    if not candidates:
        raise ValueError('No selected docking/scoring combination is executable')
    dataset_core = [{k: r.get(k) for k in ('compound_id', 'canonical_smiles', 'endpoint', 'p_value', 'value', 'units', 'relation', 'organism', 'organism_basis', 'assay', 'assay_description', 'assay_metadata', 'source_id', 'source_url', 'publication_url', 'cohort', 'raw')} for r in records]
    for row in dataset_core:
        row['measurement_id'] = row['source_id']
    review = {'schema': 2, 'identity': dataset['identity'], 'records': dataset_core,
              'record_exclusions': record_exclusions, 'ligand_subset': subset,
              'cohorts': dict(Counter(r['cohort'] for r in records)),
              'source_coverage': {k: dataset[k] for k in ('source_total', 'loaded', 'cap', 'complete', 'created_at')},
              'reference': reference, 'selection': deepcopy(selected), 'candidates': candidates, 'excluded': excluded,
              'protocol': {'scope': 'exploratory_stratified_assay_ranking', 'metrics': ['Spearman rho', 'Kendall tau', 'cognate ligand RMSD (Å)', 'PoseBusters validity', 'runtime (s)'],
                           'selection_rule': 'Display all protocols; no automatic winner or held-out performance claim',
                           'normalization': 'pX = -log10(X in mol/L); exact same-assay endpoint only',
                           'preparation': 'Meeko protein chain; any omitted incomplete terminal residues >12 Å from ligand are listed in reference; RDKit ETKDGv3 seed 42, input ligand charge/tautomer retained',
                           'biocentric_preprocessing': 'Native ensemble preprocessing with a fixed experimental crystal-ligand crop reference for every pose; 7/15 Å model-specific crops unchanged',
                           'pocket': 'Selected method and fixed parameters per chain; predicted pockets use protein only, never the reference ligand',
                           'top_pose': 'Best rescoring value among generated poses; ties retain output order',
                           'per_tool_timeout_seconds': 120, 'total_budget_seconds': None,
                           'cpu_policy': 'Per-tool reviewed thread settings', 'seed': 42, 'training_overlap': 'Unknown; public data may overlap learned-model training',
                           'construct_match': 'Assay construct/conditions and crystal construct require review; ranking is exploratory',
                           'inactive_labels': 'Not inferred', 'prep_pocket_catalog': 'Preparation selections remain UI-only; pocket methods execute in each chain',
                           'aggregation': 'Equal-weight mean of within-assay Spearman and Kendall correlations; no cross-assay pooling. Cohorts with <3 records or constant values have no ranking metric; all their records are still evaluated.'},
              'runtime': {k: os.environ.get(k) for k in ('BIOCENTRIC_VINA', 'BIOCENTRIC_GNINA', 'BIOCENTRIC_UNIDOCK', 'BIOCENTRIC_MODEL_DIR', 'BIOCENTRIC_MODEL_PYTHON')},
              'created_at': science._now()}
    from .bayesian_search import configuration
    review['search'] = configuration(candidates)
    review['protocol']['top_pose'] += '; execution is ligand-first within each Bayesian search round'
    review['protocol']['selection_rule'] = 'Automatically fix the best available evaluated chain; failed validation gates remain explicit; no held-out performance claim'
    review['runtime_fingerprint'] = runtime_fingerprint(candidates)
    review['sha256'] = digest(review)
    return review


def validate_review(review):
    if digest({k: v for k, v in review.items() if k != 'sha256'}) != review.get('sha256'):
        raise ValueError('Evaluation plan changed after review')
    if review.get('schema', 1) < 2 and len({r['cohort'] for r in review['records']}) != 1:
        raise ValueError('Different assay cohorts cannot share a benchmark')
    if review.get('schema', 1) >= 2 and len({measurement_key(r) for r in review['records']}) != len(review['records']):
        raise ValueError('Duplicate measurement identities in evaluation set')
    if any(not math.isfinite(r['p_value']) for r in review['records']):
        raise ValueError('Non-finite experimental labels')
    if review.get('search'):
        from .bayesian_search import configuration
        config = configuration(review['candidates'], review['search'].get('method'))
        if review['search'] != config:
            raise ValueError('Unsupported Bayesian search configuration')


def start_run(root, review):
    validate_review(review)
    if runtime_fingerprint(review['candidates']) != review['runtime_fingerprint']:
        raise ValueError('Runtime changed after review. Prepare a new evaluation plan.')
    folder = Path(root).resolve() / 'evaluation-runs' / uuid.uuid4().hex
    folder.mkdir(parents=True)
    with (folder / 'plan.json').open('x') as f:
        json.dump(review, f, indent=2, allow_nan=False)
    science._write_json(folder / 'status.json', {'status': 'queued', 'rows': [], 'events': [], 'stage': 'queued', 'completed': 0,
                                               'total': (review.get('search') or {}).get('max_trials', len(review['candidates'])), 'progress_unit': 'combinations', 'plan_sha256': review['sha256']})
    with (folder / 'worker.log').open('w') as log:
        child = subprocess.Popen([sys.executable, '-m', 'biocentric.discovery_benchmark', 'worker', str(folder)],
                                 cwd=Path(__file__).resolve().parents[1], stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                 start_new_session=True, env=science._tool_environment())
    science._write_json(folder / 'process.json', {'pid': child.pid, 'created_at': science._now()})
    return str(folder)


def read_status(folder):
    folder = Path(folder)
    path = folder / 'status.json'
    state = json.loads(path.read_text())
    process_file = folder / 'process.json'
    if (state['status'] in {'queued', 'running'} and process_file.is_file() and Path('/proc').is_dir()
            and time.time() - path.stat().st_mtime > 10):
        pid = json.loads(process_file.read_text())['pid']
        try:
            command = (Path('/proc') / str(pid) / 'cmdline').read_bytes()
        except OSError:
            command = b''
        if b'biocentric.discovery_benchmark' not in command or str(folder.resolve()).encode() not in command:
            state.update(status='interrupted', stage='stopped', error='Evaluation worker is no longer active. Completed results are preserved; unfinished pairs have no score.')
    return state


def cancel_run(folder):
    (Path(folder) / 'cancel.request').touch()


def summary_metrics(review, rows):
    from .pipeline_decisions import stratified_statistics
    summaries = []
    for pair in review['candidates']:
        good = [r for r in rows if r['combination'] == pair['id'] and r.get('kind') == 'measured' and r.get('status') == 'completed']
        summary = {'Combination': chain_label(pair), 'Measurements': len(good), 'Expected': len(review['records']), 'Spearman ρ': None, 'Kendall τ': None}
        if len(good) == len(review['records']) and len(good) >= 3:
            sign = 1 if SCORES[pair['scoring']['tool']][2] else -1
            stats = stratified_statistics(good, sign, intervals=False)
            summary.update({'Spearman ρ': stats['rho'], 'Kendall τ': stats['tau'], 'Assays ranked': stats['ranked_cohorts']})
        summaries.append(summary)
    return summaries


def _write_sdf(path, mol):
    from rdkit import Chem
    path.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    return str(path)


def _pdbqt(mol, path):
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    setups = MoleculePreparation(charge_model='gasteiger').prepare(Chem.AddHs(mol, addCoords=True))
    if len(setups) != 1:
        raise ValueError('Ambiguous ligand preparation')
    text, success, error = PDBQTWriterLegacy.write_string(setups[0])
    if not success:
        raise ValueError(error)
    path.write_text(text)
    return str(path)


def _dock(pair, mol, folder, receptor, pdbqt, box, cancelled, emit):
    from .engines import run_engine
    from .engine_workflow import _native_molecules
    folder.mkdir(parents=True)
    engine = pair['docking']['tool']
    sdf = _write_sdf(folder / 'input.sdf', mol)
    par = {**pair['docking']['parameters'], 'box': {k: box[k] for k in ('center', 'size')}, 'timeout_seconds': 120}
    if engine == 'gnina':
        par.setdefault('device', 'cpu')  # Legacy plans retain their reviewed CPU default.
        inputs = {'receptor': str(receptor), 'ligand': sdf, 'preparation_reviewed': True}
    else:
        ligand = _pdbqt(mol, folder / 'input.pdbqt')
        inputs = {'receptor': str(pdbqt), 'preparation_reviewed': True, **({'ligands': [ligand]} if engine == 'unidock' else {'ligand': ligand})}
    native = run_engine(engine, inputs, par, folder / 'native', emit=emit, cancelled=cancelled)
    if native['status'] != 'completed':
        raise ValueError('; '.join(native.get('blockers', [])) or native.get('error', {}).get('message') or 'Docking failed')
    poses = _native_molecules(engine, native)
    if not poses or any(p is None or science._canonical(p) != science._canonical(mol) for p in poses):
        raise ValueError('Docked chemistry differs from the measured compound')
    return poses


def _score(pair, poses, folder, receptor, pdbqt, cancelled, emit, box=None, archive=None):
    from .engines import run_engine, _executable
    from rdkit import Chem
    score = pair['scoring']['tool']
    folder.mkdir(parents=True, exist_ok=True)
    if score == 'biocentric_affinity':
        from .discovery_affinity import score_groups
        group = score_groups({'poses': poses}, receptor, folder, pair['scoring']['parameters'], cancelled, emit)['poses']
        if group['errors']:
            raise ValueError('; '.join(group['errors']))
        return group['scores']
    if score.startswith('gnina'):
        source = folder / 'poses.sdf'
        with Chem.SDWriter(str(source)) as writer:
            for mol in poses:
                writer.write(mol)
        device = pair['scoring']['parameters'].get('device', 'cpu')
        parameters = {'cpu': pair['scoring']['parameters'].get('cpu', 1), 'device': device, 'timeout_seconds': 120}
        def rescore(out):
            out.mkdir(parents=True, exist_ok=True)
            result = run_engine('gnina_affinity', {'receptor': str(receptor), 'ligand': str(source), 'preparation_reviewed': True},
                                parameters, out / 'native', emit=emit, cancelled=cancelled)
            if result['status'] != 'completed':
                raise ValueError('; '.join(result.get('blockers', [])) or result.get('error', {}).get('message') or 'Rescoring failed')
            output = list(Chem.SDMolSupplier(result['expected_outputs'][0], removeHs=False))
            if len(output) != len(poses):
                raise ValueError('Scoring output count differs from docked poses')
            for before, after in zip(poses, output):
                check = science._symmetry_rmsd(before, after)
                if check.get('rmsd_angstrom') is None or check['rmsd_angstrom'] > .02:
                    raise ValueError('Scoring changed pose geometry or identity')
            metrics = {key: [next(m['value'] for m in r['metrics'] if m['name'] == key) for r in result['results']]
                       for key in ('CNNscore', 'CNNaffinity')}
            if any(len(values) != len(poses) or any(not math.isfinite(v) for v in values) for values in metrics.values()):
                raise ValueError('Missing or nonfinite GNINA scores')
            science._write_json(out / 'scores.json', metrics)
        destination = folder / 'fixed-pose-scores'
        if archive:
            archive.materialize('gnina_scoring', {'tool': 'gnina', 'receptor': science.sha256(receptor),
                'poses': science.sha256(source), 'parameters': parameters}, destination, rescore)
        else:
            rescore(destination)
        key = 'CNNscore' if score == 'gnina_cnn' else 'CNNaffinity'
        return json.loads((destination / 'scores.json').read_text())[key]
    values = []
    for i, mol in enumerate(poses):
        ligand = _pdbqt(mol, folder / f'pose-{i}.pdbqt')
        log = folder / f'score-{i}.log'
        grid = box
        if grid is None:
            xyz = mol.GetConformer().GetPositions()
            grid = {'center': ((xyz.min(axis=0) + xyz.max(axis=0)) / 2).tolist(), 'size': (xyz.max(axis=0) - xyz.min(axis=0) + 8).tolist()}
        args = []
        for dimension, center, size in zip('xyz', grid['center'], grid['size']):
            args += [f'--center_{dimension}', str(center), f'--size_{dimension}', str(size)]
        science.run_process([_executable('vina'), '--receptor', str(pdbqt), '--ligand', ligand, '--score_only', '--cpu', str(pair['scoring']['parameters'].get('cpu', 1)), *args], folder, log, 120, cancelled)
        match = re.search(r'(?:Affinity|Estimated Free Energy of Binding)\s*:\s*([-+\d.eE]+)', log.read_text())
        if not match or not math.isfinite(float(match[1])):
            raise ValueError('No finite Vina score returned')
        values.append(float(match[1]))
    return values


def _prepare_docked_poses(pair, record, work, receptor, pdbqt, box, cache, cancelled, emit):
    from rdkit import Chem
    method = pair.get('pocket', pockets.REFERENCE)
    if isinstance(box, str):
        raise ValueError('Pocket identification failed: ' + box)
    cache_key = digest([method, box, pair['docking'], record['canonical_smiles']])
    archive = getattr(cache, 'archive', None)
    if cache_key not in cache:
        work.mkdir(parents=True, exist_ok=True)
        ligand_key = ('ligand', record['canonical_smiles'])
        if ligand_key not in cache:
            def prepare(out):
                out.mkdir(parents=True, exist_ok=True)
                inp = out / 'input.smi'
                inp.write_text(record['canonical_smiles'] + ' ligand\n')
                lib = science.prepare_library(inp, out / 'ligand', seed=42, cancelled=cancelled)
                if not lib['compounds']:
                    raise ValueError(str(lib['exclusions']))
                mol = next(iter(Chem.SDMolSupplier(lib['prepared_sdf'], removeHs=False)))
                _write_sdf(out / 'prepared.sdf', mol)
            destination = work / 'prepared-ligand'
            if archive:
                archive.materialize('ligand', {'smiles': record['canonical_smiles'], 'seed': 42}, destination, prepare)
            else:
                prepare(destination)
            cache[ligand_key] = next(iter(Chem.SDMolSupplier(str(destination / 'prepared.sdf'), removeHs=False)))
        mol = cache[ligand_key]
        if archive:
            def dock(out):
                out.mkdir(parents=True, exist_ok=True)
                poses = _dock(pair, mol, out / 'dock', receptor, pdbqt, box, cancelled, emit)
                with Chem.SDWriter(str(out / 'poses.sdf')) as writer:
                    for pose in poses:
                        writer.write(pose)
            destination = work / 'docked-poses'
            archive.materialize('docking', {'tool': pair['docking']['tool'], 'configuration': pair['docking'],
                'receptor': science.sha256(receptor), 'prepared_receptor': science.sha256(pdbqt),
                'pocket': method, 'box': box, 'ligand': Chem.MolToMolBlock(mol)}, destination, dock)
            poses = list(Chem.SDMolSupplier(str(destination / 'poses.sdf'), removeHs=False))
            import numpy as np
            if not poses or any(p is None or science._canonical(p) != science._canonical(mol)
                                or not np.isfinite(p.GetConformer().GetPositions()).all() for p in poses):
                raise ValueError('Archived docking poses failed integrity validation')
            cache[cache_key] = poses
        else:
            cache[cache_key] = _dock(pair, mol, work / 'dock', receptor, pdbqt, box, cancelled, emit)
    else:
        emit({'tool': pair['docking']['tool'], 'stage': 'archive', 'message': 'Reused docking poses across scoring chains'})
    return cache[cache_key]


def _affinity_batch(pair, records, start, folder, receptor, pdbqt, box, cache,
                    reference, include_reference, state, save, cancelled, emit):
    from rdkit import Chem
    from .discovery_affinity import score_groups
    groups, errors, durations = {}, {}, {}
    for index, record in enumerate(records, start):
        science._cancel(cancelled)
        state.update(active_compound=record['compound_id'], stage='docking')
        save()
        started = time.monotonic()
        try:
            groups[str(index)] = _prepare_docked_poses(pair, record, folder / 'calculations' / pair['id'] / str(index),
                                                     receptor, pdbqt, box, cache, cancelled, emit)
        except Exception as exc:
            science._cancel(cancelled)
            errors[str(index)] = {'errors': [str(exc)]}
        durations[str(index)] = time.monotonic() - started
    if include_reference and not isinstance(box, str):
        groups['crystal_reference'] = [Chem.AddHs(reference, addCoords=True)]
    if groups:
        state.update(stage='scoring', active_compound=f'Batch {start + 1}–{start + len(records)}')
        save()
        scoring_started = time.monotonic()
        try:
            results = score_groups(groups, receptor, folder / 'affinity-batches' / pair['id'] / str(start),
                                   pair['scoring']['parameters'], cancelled, emit)
        except Exception as exc:
            science._cancel(cancelled)
            results = {key: {'errors': [str(exc)]} for key in groups}
        errors.update(results)
        # Attribute shared batch inference time equally; preserve native batch
        # duration/resources in result.json for direct runtime inspection.
        share = (time.monotonic() - scoring_started) / len(records)
        durations = {key: value + share for key, value in durations.items()}
    for key, duration in durations.items():
        errors[key]['runtime_seconds'] = duration
    return errors


def worker(folder):
    from rdkit import Chem
    folder = Path(folder).resolve()
    review = json.loads((folder / 'plan.json').read_text())
    validate_review(review)
    if review.get('search'):
        from .discovery_search_runner import worker as search_worker
        return search_worker(folder, review)
    state = read_status(folder)
    started = time.monotonic()
    state.update(status='running', started_at=science._now(), stage='preparation')
    state.update(active_combination=review['candidates'][0]['id'], active_compound='Receptor')
    def save():
        state['updated_at'] = science._now()
        state['elapsed_seconds'] = time.monotonic() - started
        if review.get('schema', 1) >= 2:
            progress = combination_progress(review, state)
            state.update(completed=progress['completed'], total=progress['total'], progress_unit='combinations')
        science._write_json(folder / 'status.json', state)
    def cancelled():
        budget = review['protocol'].get('total_budget_seconds')
        return (folder / 'cancel.request').exists() or (budget is not None and time.monotonic() - started > budget)
    def emit(event):
        if event.get('stage') == 'archive':
            counts = state.setdefault('reuse', {'hits': 0, 'stored': 0})
            counts['hits' if event.get('message', '').startswith('Reused') else 'stored'] += 1
        state['events'] = (state['events'] + [{k: event.get(k) for k in ('tool', 'stage', 'message')}])[-40:]
        save()
    try:
        if runtime_fingerprint(review['candidates']) != review['runtime_fingerprint']:
            raise ValueError('Runtime changed after approval')
        save()
        receptor = folder / 'receptor.pdb'
        reference_path = folder / 'reference.sdf'
        receptor.write_text(review['reference']['receptor_pdb'])
        reference_path.write_text(review['reference']['reference_sdf'])
        prepared = science.prepare_receptor(receptor, folder / 'receptor-preparation', 'meeko', reviewed=True, cancelled=cancelled)
        pdbqt = Path(prepared['output'])
        reference = next(iter(Chem.SDMolSupplier(str(reference_path), removeHs=False)))
        tasks = [{'compound_id': 'crystal:' + review['reference']['component'], 'canonical_smiles': review['reference']['canonical_smiles'], 'kind': 'reference', 'p_value': None, 'endpoint': None},
                 *[{**r, 'kind': 'measured'} for r in review['records']]]
        from .evaluation_archive import Archive, StageCache, identify_pocket
        archive = (Archive(folder.parent.parent / 'evaluation-archive', review['runtime_fingerprint'], cancelled, emit)
                   if review['runtime_fingerprint'].get('files') else None)
        cache, reference_scores, pocket_cache = StageCache(archive), {}, {}
        for pair in review['candidates']:
            method = pair.get('pocket', pockets.REFERENCE)
            pocket_key = digest(method)
            state.update(active_combination=pair['id'], active_compound='Pocket', stage='pocket')
            save()
            if pocket_key not in pocket_cache:
                try:
                    pocket_cache[pocket_key] = identify_pocket(method, receptor, reference_path, folder / 'pockets' / pocket_key[:16], cancelled, archive)
                except Exception as exc:
                    if cancelled():
                        raise
                    pocket_cache[pocket_key] = str(exc)
            box = pocket_cache[pocket_key]
            affinity = pair['scoring']['tool'] == 'biocentric_affinity'
            batch_results, crystal_prediction = {}, None
            if affinity:
                from .discovery_affinity import Parameters
                batch_size = Parameters.model_validate(pair['scoring']['parameters']).batch_size
            for index, record in enumerate(tasks):
                science._cancel(cancelled)
                if affinity and index % batch_size == 0:
                    batch_results = _affinity_batch(pair, tasks[index:index + batch_size], index, folder,
                        receptor, pdbqt, box, cache, reference, index == 0, state, save, cancelled, emit)
                    if index == 0:
                        crystal_prediction = batch_results.get('crystal_reference')
                state.update(active_combination=pair['id'], active_compound=record['compound_id'], stage='docking')
                save()
                row = {'combination': pair['id'], 'docking': pair['docking']['tool'], 'scoring': pair['scoring']['tool'],
                       'compound_id': record['compound_id'], 'measurement_id': measurement_key(record), 'cohort': record.get('cohort'), 'pocket': method['tool'], 'kind': record['kind'], 'p_value': record['p_value'],
                       'endpoint': record['endpoint'], 'score': None, 'score_unit': SCORES[pair['scoring']['tool']][1],
                       'rmsd_angstrom': None, 'physical_pass': None, 'status': 'failed'}
                t = time.monotonic()
                try:
                    work = folder / 'calculations' / pair['id'] / str(index)
                    work.mkdir(parents=True, exist_ok=True)
                    if isinstance(box, str):
                        raise ValueError('Pocket identification failed: ' + box)
                    row['pocket_box'] = box
                    prediction = batch_results.get(str(index)) if affinity else None
                    if affinity and (not prediction or prediction['errors']):
                        raise ValueError('; '.join(prediction['errors']) if prediction else 'Missing BioCentric pose group')
                    poses = _prepare_docked_poses(pair, record, work, receptor, pdbqt, box, cache, cancelled, emit)
                    state['stage'] = 'scoring'
                    save()
                    scores = prediction['scores'] if affinity else _score(pair, poses, work / 'score', receptor, pdbqt, cancelled, emit, box)
                    if len(scores) != len(poses) or any(not math.isfinite(v) for v in scores):
                        raise ValueError('Missing or nonfinite scores')
                    best = (max if SCORES[pair['scoring']['tool']][2] else min)(range(len(scores)), key=lambda i: scores[i])
                    pose = poses[best]
                    pose_path = Path(_write_sdf(work / 'best-pose.sdf', pose))
                    state['stage'] = 'validation'
                    save()
                    validity = science.validate_pose(pose_path, receptor, work / 'validation', 120, cancelled)
                    same = science._canonical(pose) == science._canonical(reference)
                    score_key = digest([method, box, pair['scoring']])
                    if same and score_key not in reference_scores:
                        if affinity:
                            if not crystal_prediction or crystal_prediction['errors']:
                                raise ValueError('BioCentric crystal reference scoring failed: ' + '; '.join((crystal_prediction or {}).get('errors', ['Missing reference'])))
                            reference_scores[score_key] = crystal_prediction['scores'][0]
                        else:
                            reference_scores[score_key] = _score(pair, [Chem.AddHs(reference, addCoords=True)], work / 'reference-score', receptor, pdbqt, cancelled, emit, box)[0]
                    if affinity:
                        row['model_prediction'] = {
                            'device': prediction['device'], 'provenance': prediction['provenance'],
                            'best_pose': prediction['predictions'][best],
                            'artifact': str(Path(prediction['artifact']).relative_to(folder)),
                            'artifact_sha256': science.sha256(Path(prediction['artifact'])),
                            'ensemble_spread_calibrated': False, 'experimental_measurement': False}
                    row.update(status='completed', score=scores[best], pose_path=str(pose_path.relative_to(folder)),
                               pose_sha256=science.sha256(pose_path), physical_pass=validity.get('passed'),
                               physical_status=validity.get('status'),
                               rmsd_angstrom=science._symmetry_rmsd(pose, reference).get('rmsd_angstrom') if same else None,
                               reference_comparable=same, reference_score=reference_scores.get(score_key) if same else None,
                               top_pose_index=best + 1, all_scores=scores)
                except Exception as exc:
                    if cancelled():
                        raise science.ScienceError('cancelled', 'Run stopped or reached its configured time budget') from exc
                    row['error'] = str(exc)
                row['runtime_seconds'] = time.monotonic() - t + (batch_results.get(str(index), {}).get('runtime_seconds', 0) if affinity else 0)
                state['rows'].append(row)
                if review.get('schema', 1) < 2:
                    state['completed'] += 1
                if index == len(tasks) - 1:
                    state['metrics'] = summary_metrics(review, state['rows'])
                save()
        failed = sum(r['status'] != 'completed' for r in state['rows'])
        state.update(status='completed' if not failed else 'partial' if failed < len(state['rows']) else 'failed', stage='finished')
    except Exception as exc:
        state.update(status='cancelled' if cancelled() else 'failed', error=str(exc), stage='stopped')
    finally:
        # Every unexecuted pair stays accounted for, including cancellation/fatal preparation failure.
        accounted = {(r['combination'], measurement_key(r)) for r in state['rows']}
        for pair in review['candidates']:
            for record in [{'compound_id': 'crystal:' + review['reference']['component'], 'kind': 'reference', 'p_value': None, 'endpoint': None}, *[{**r, 'kind': 'measured'} for r in review['records']]]:
                if (pair['id'], measurement_key(record)) not in accounted:
                    state['rows'].append({'combination': pair['id'], 'docking': pair['docking']['tool'], 'scoring': pair['scoring']['tool'], 'compound_id': record['compound_id'], 'measurement_id': measurement_key(record), 'cohort': record.get('cohort'), 'pocket': pair.get('pocket', pockets.REFERENCE)['tool'], 'kind': record['kind'], 'p_value': record['p_value'], 'endpoint': record['endpoint'], 'status': 'not_run', 'score': None, 'error': state.get('error', 'Not run')})
        state['metrics'] = summary_metrics(review, state['rows'])
        save()


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == 'worker':
        worker(sys.argv[2])
