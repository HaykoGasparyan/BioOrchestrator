"""Public target evidence for the discovery workspace; never runs docking jobs."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import uuid

import httpx

from .retrieval import RetrievalClient, RetrievalError, ACCESSION
from .viewer import VIEWER_SHA256, VIEWER_VERSION, script_json

PAGE_SIZE = 24
ENTRY_QUERY = '''query($ids:[String!]!){entries(entry_ids:$ids){
 rcsb_id struct{title} exptl{method}
 rcsb_entry_info{resolution_combined deposited_model_count}
 rcsb_primary_citation{title year pdbx_database_id_PubMed pdbx_database_id_DOI}
 polymer_entities{rcsb_polymer_entity{pdbx_description pdbx_mutation}
 rcsb_polymer_entity_container_identifiers{auth_asym_ids uniprot_ids}}
 nonpolymer_entities{nonpolymer_comp{chem_comp{id name formula_weight}}}
}}'''


def client(root):
    c = RetrievalClient(Path(root) / 'discovery-cache')
    c.timeout, c.max_attempts, c.min_interval = 12.0, 1, .25
    return c


def resolve_identity(target, organism, root):
    if re.fullmatch(r'[0-9][A-Za-z0-9]{3}', target.strip()):
        pdb_id = target.strip().upper()
        c = client(root)
        entry = c._request('RCSB PDB', 'https://data.rcsb.org/graphql', payload={'query': ENTRY_QUERY, 'variables': {'ids': [pdb_id]}})
        entries = entry['data'].get('data', {}).get('entries') or []
        accessions = sorted({accession for row in entries if row for entity in row.get('polymer_entities') or []
                             for accession in entity.get('rcsb_polymer_entity_container_identifiers', {}).get('uniprot_ids') or []})
        candidates = []
        for accession in accessions:
            resolved = c.resolve_target(accession, organism)
            candidates.extend({**identity, 'requested_pdb': pdb_id} for identity in resolved['candidates'])
        return {'status': 'resolved' if len(candidates) == 1 else 'ambiguous' if candidates else 'not_found',
                'candidates': candidates, 'requested_pdb': pdb_id, 'limitations': []}
    result = client(root).resolve_target(target, organism)
    exact = [r for r in result['candidates'] if r.get('reviewed') and target.casefold() in {str(g).casefold() for g in r.get('gene_names', [])}]
    if len(exact) == 1:
        result.update(status='resolved', candidates=exact)
    return result


def _accession(identity):
    accession = identity.get('accession', '')
    if not ACCESSION.fullmatch(accession) or '-' in accession:
        raise RetrievalError('A canonical UniProt accession is required.')
    return accession


def structure_page(identity, root, page=0, entry_ids=None):
    accession = _accession(identity)
    if not isinstance(page, int) or page < 0:
        raise ValueError('Invalid page')
    c = client(root)
    query = {'query': {'type': 'terminal', 'service': 'text', 'parameters': {
        'attribute': 'rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession',
        'operator': 'exact_match', 'value': accession}}, 'return_type': 'entry',
        'request_options': {'paginate': {'start': page * PAGE_SIZE, 'rows': PAGE_SIZE},
                            'results_content_type': ['experimental'], 'sort': [{'sort_by': 'rcsb_id', 'direction': 'asc'}]}}
    if entry_ids is not None:
        ids = list(dict.fromkeys(str(value).upper() for value in entry_ids))
        if any(not re.fullmatch(r'[0-9][A-Z0-9]{3}', value) for value in ids):
            raise ValueError('Invalid explicit PDB identifier')
        hits = {'data': {'total_count': len(ids)}, 'retrieved_at': None}
    else:
        hits = c._request('RCSB PDB', 'https://search.rcsb.org/rcsbsearch/v2/query', payload=query)
        ids = [row['identifier'] for row in hits['data'].get('result_set', [])
               if re.fullmatch(r'[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]{8}', row.get('identifier', ''))]
    if not ids:
        return {'rows': [], 'total': hits['data'].get('total_count', 0), 'page': page, 'retrieved_at': hits['retrieved_at']}
    entries = c._request('RCSB PDB', 'https://data.rcsb.org/graphql', payload={'query': ENTRY_QUERY, 'variables': {'ids': ids}})
    rows = []
    for row in entries['data'].get('data', {}).get('entries', []) or []:
        if not row or row.get('rcsb_id') not in ids:
            continue
        target_entities = [p for p in row.get('polymer_entities') or []
                           if accession in (p.get('rcsb_polymer_entity_container_identifiers', {}).get('uniprot_ids') or [])]
        # Verify returned metadata maps to the requested protein, not just a text match.
        if not target_entities:
            continue
        compounds = [p.get('nonpolymer_comp', {}).get('chem_comp', {}) for p in row.get('nonpolymer_entities') or []]
        compounds = [p for p in compounds if p.get('id') and p['id'] not in {'HOH', 'DOD'}]
        info = row.get('rcsb_entry_info') or {}
        rows.append({'id': row['rcsb_id'], 'title': (row.get('struct') or {}).get('title', ''),
                     'methods': [e['method'] for e in row.get('exptl') or []],
                     'resolution': (info.get('resolution_combined') or [None])[0],
                     'models': info.get('deposited_model_count') or 1,
                     'compounds': compounds,
                     'partners': [{'name': p.get('rcsb_polymer_entity', {}).get('pdbx_description') or 'Other polymer', 'chains': p.get('rcsb_polymer_entity_container_identifiers', {}).get('auth_asym_ids') or []} for p in row.get('polymer_entities') or [] if p not in target_entities],
                     'target_chains': [ch for p in target_entities for ch in p['rcsb_polymer_entity_container_identifiers'].get('auth_asym_ids') or []],
                     'mutations': [p['rcsb_polymer_entity'].get('pdbx_mutation') for p in target_entities if p.get('rcsb_polymer_entity', {}).get('pdbx_mutation')],
                     'citation': row.get('rcsb_primary_citation') or {},
                     'url': f"https://www.rcsb.org/structure/{row['rcsb_id']}", 'retrieved_at': entries['retrieved_at']})
    rows.sort(key=lambda r: (not any((p.get('formula_weight') or 0) > 150 for p in r['compounds']), r['id']))
    return {'rows': rows, 'total': hits['data'].get('total_count', 0), 'page': page,
            'retrieved_at': entries['retrieved_at'], 'matched_on_page': len(ids)}


def binding_page(identity, root, page=0):
    accession = _accession(identity)
    if not isinstance(page, int) or page < 0:
        raise ValueError('Invalid page')
    c = client(root)
    targets = c._request('ChEMBL', 'https://www.ebi.ac.uk/chembl/api/data/target.json',
                        {'target_components__accession': accession, 'target_type': 'SINGLE PROTEIN', 'limit': 20})
    target_ids = [t['target_chembl_id'] for t in targets['data'].get('targets', [])
                  if str(t.get('tax_id')) == str(identity.get('taxon_id')) and t.get('target_type') == 'SINGLE PROTEIN'
                  and any(p.get('accession') == accession for p in t.get('target_components', []))]
    if not target_ids:
        return {'rows': [], 'total': 0, 'page': page, 'retrieved_at': targets['retrieved_at']}
    data = c._request('ChEMBL', 'https://www.ebi.ac.uk/chembl/api/data/activity.json',
                      {'target_chembl_id__in': ','.join(target_ids), 'assay_type': 'B', 'standard_type__in': 'Kd,Ki,IC50',
                       'limit': 50, 'offset': 50 * page, 'order_by': 'activity_id'})
    records = [r for r in data['data'].get('activities', []) if r.get('target_chembl_id') in target_ids and r.get('assay_type') == 'B']
    publications = sorted({r['document_chembl_id'] for r in records if re.fullmatch(r'CHEMBL[0-9]+', r.get('document_chembl_id') or '')})
    docs = {}
    if publications:
        try:
            doc_data = c._request('ChEMBL', 'https://www.ebi.ac.uk/chembl/api/data/document.json',
                                 {'document_chembl_id__in': ','.join(publications), 'limit': 50})
            docs = {r['document_chembl_id']: r for r in doc_data['data'].get('documents', [])}
        except RetrievalError:
            pass  # The original document ID is still retained and linked below.
    rows = []
    for r in records:
        d = docs.get(r.get('document_chembl_id'), {})
        publication_url = f"https://pubmed.ncbi.nlm.nih.gov/{d['pubmed_id']}/" if str(d.get('pubmed_id', '')).isdigit() else f"https://www.ebi.ac.uk/chembl/explore/document/{r.get('document_chembl_id', '')}"
        rows.append({'compound': r.get('molecule_pref_name') or r.get('molecule_chembl_id'),
                     'compound_id': r.get('molecule_chembl_id'), 'smiles': r.get('canonical_smiles'),
                     'endpoint': r.get('standard_type'), 'relation': r.get('standard_relation'),
                     'value': r.get('standard_value'), 'units': r.get('standard_units'),
                     'assay': r.get('assay_chembl_id'), 'assay_description': r.get('assay_description'),
                     'organism': r.get('assay_organism'), 'target_id': r.get('target_chembl_id'),
                     'document': r.get('document_chembl_id'), 'publication': d.get('title') or r.get('document_chembl_id'),
                     'publication_url': publication_url, 'year': r.get('document_year'),
                     'validity': r.get('data_validity_comment'), 'comment': r.get('activity_comment'),
                     'source_url': f"https://www.ebi.ac.uk/chembl/explore/activities/ACTIVITY_ID:{r.get('activity_id')}",
                     'retrieved_at': data['retrieved_at'], 'raw': r})
    return {'rows': rows, 'total': data['data'].get('page_meta', {}).get('total_count', len(rows)), 'page': page, 'retrieved_at': data['retrieved_at']}


def literature(identity, root):
    c = client(root)
    result = {'evidence': [], 'limitations': []}
    c._literature(identity, '', result)
    return result


def overview(identity, root):
    results = {}
    with ThreadPoolExecutor(max_workers=3) as pool:
        pending = {name: pool.submit(fn, identity, root) for name, fn in
                   [('structures', structure_page), ('binding', binding_page), ('literature', literature)]}
        for name, future in pending.items():
            try:
                results[name] = future.result()
            except (RetrievalError, KeyError, TypeError, ValueError) as exc:
                results[name] = {'rows': [], 'total': None, 'error': str(exc)}
    if identity.get('requested_pdb'):
        requested = structure_page(identity, root, entry_ids=[identity['requested_pdb']])
        page = results['structures']
        page['rows'] = requested['rows'] + [row for row in page.get('rows', []) if row['id'] != identity['requested_pdb']]
    return results


def coordinates(pdb_id, root, format_name='cif'):
    if not re.fullmatch(r'[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]{8}', pdb_id):
        raise ValueError('Invalid PDB identifier')
    if format_name not in {'cif', 'pdb'}:
        raise ValueError('Unsupported coordinate format')
    folder = Path(root) / 'discovery-cache' / 'coordinates'
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / (pdb_id + ('.pdb' if format_name == 'pdb' else '') + '.json')
    if path.exists():
        try:
            cached = json.loads(path.read_text())
            if cached.get('id') == pdb_id and hashlib.sha256(cached['text'].encode()).hexdigest() == cached.get('sha256'):
                cached['format'] = format_name
                return cached
        except (ValueError, KeyError):
            pass
    url = f'https://files.rcsb.org/download/{pdb_id}.{format_name}'
    with httpx.Client(timeout=20, follow_redirects=False) as session:
        with session.stream('GET', url) as response:
            response.raise_for_status()
            body = bytearray()
            for chunk in response.iter_bytes():
                body.extend(chunk)
                if len(body) > 12_000_000:
                    raise RetrievalError('This structure exceeds the 12 MB interactive preview limit. Use its PDB link.')
    text = body.decode('utf-8')
    if (format_name == 'cif' and '_atom_site.' not in text) or (format_name == 'pdb' and not any(line.startswith(('ATOM  ', 'HETATM')) for line in text.splitlines())):
        raise RetrievalError('The downloaded record contains no experimental atom coordinates.')
    result = {'id': pdb_id, 'format': format_name, 'text': text, 'sha256': hashlib.sha256(text.encode()).hexdigest(),
              'url': url, 'retrieved_at': datetime.now(timezone.utc).isoformat()}
    temp = path.with_suffix('.tmp.' + uuid.uuid4().hex)
    temp.write_text(json.dumps(result))
    temp.replace(path)
    return result


def experimental_viewer(record, geometry, ligand=''):
    if record['id'] != geometry['id'] or hashlib.sha256(geometry['text'].encode()).hexdigest() != geometry['sha256']:
        raise ValueError('Experimental structure and coordinates do not match.')
    library = (Path(__file__).parent / 'assets' / f'3Dmol-{VIEWER_VERSION}-min.js').read_bytes()
    if hashlib.sha256(library).hexdigest() != VIEWER_SHA256:
        raise ValueError('3D viewer integrity check failed.')
    data = {'text': geometry['text'], 'format': geometry.get('format', 'cif'), 'ligand': ligand, 'chains': record['target_chains'], 'expected_models': record.get('models', 1)}
    return '''<!doctype html><html><head><meta charset="utf-8"><meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src data: blob:; connect-src 'none'; worker-src 'none'; object-src 'none'; base-uri 'none'; form-action 'none'">
<style>body{margin:0;font:12px system-ui;color:#536681;background:#f5f7fc}#mol{height:370px;position:relative}.bar{display:flex;gap:7px;align-items:center;padding:9px 12px}button,select{border:1px solid #dce2f2;border-radius:7px;padding:7px 10px;background:white;color:#5267a4;cursor:pointer}.hint{margin-left:auto;color:#8794ad}#error{padding:12px;color:#b84c59}</style></head><body>
<div class="bar"><button id="fit">Structure</button><button id="site">Ligand</button><button id="spin">Rotate</button><select id="frame" aria-label="Deposited model"></select><span class="hint">Drag · scroll to zoom</span></div><div id="mol"></div><div id="error"></div>
<script>''' + library.decode().replace('</script', '<\\/script') + '''</script><script type="application/json" id="data">''' + script_json(data) + '''</script><script>
(async function(){try{const d=JSON.parse(document.getElementById('data').textContent),v=$3Dmol.createViewer('mol',{backgroundColor:'#f5f7fc',antialias:true});
v.addModelsAsFrames(d.text,d.format);if(!v.getModel()||!v.getModel().selectedAtoms({}).length)throw Error('No atoms could be displayed.');
const ligand=d.ligand.startsWith('chain:')?{chain:d.ligand.slice(6)}:d.ligand?{resn:d.ligand}:{hetflag:true,not:{resn:['HOH','DOD']}};
function style(){v.setStyle({},{cartoon:{color:'#bbc6df',opacity:.55}});if(d.chains.length)v.setStyle({chain:d.chains},{cartoon:{color:'spectrum',opacity:.95}});v.setStyle({hetflag:true,not:{resn:['HOH','DOD']}},{stick:{radius:.15,colorscheme:'cyanCarbon'}});if(d.ligand)v.setStyle(ligand,{stick:{radius:.23,colorscheme:'orangeCarbon'}});v.render();}
style();v.zoomTo();v.render();let spin=false;document.getElementById('fit').onclick=()=>{v.zoomTo();v.render()};document.getElementById('site').disabled=!v.getModel().selectedAtoms(ligand).length;document.getElementById('site').onclick=()=>{v.zoomTo(ligand);v.zoom(.75);v.render()};document.getElementById('spin').onclick=()=>{spin=!spin;v.spin(spin?'y':false)};
const select=document.getElementById('frame'),n=v.getNumFrames();if(n!==d.expected_models)document.getElementById('error').textContent='Showing '+n+' of '+d.expected_models+' deposited models. Download the coordinates for full inspection.';select.hidden=n<2;for(let i=0;i<n;i++){const o=document.createElement('option');o.value=i;o.textContent='Model '+(i+1)+' / '+n;select.appendChild(o)}select.onchange=async()=>{await v.setFrame(Number(select.value));style()};window.addEventListener('resize',()=>{v.resize();v.render()});
}catch(e){document.getElementById('error').textContent='3D unavailable: '+e.message}})();</script></body></html>'''
