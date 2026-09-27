"""Adapter inventory separates implementation, installed health and validation."""
import importlib.metadata
import importlib.util
import os
import shutil
from pathlib import Path
from .models import ToolSpec


def adapters() -> list[dict]:
    source_data=[
        ("UniProt","Target identity and sequence","https://www.uniprot.org/help/api_queries","REST v1","protein query → identity"),
        ("Europe PMC","Abstract evidence","https://europepmc.org/RestfulWebService","REST v6.9","query → abstract records"),
        ("Open Targets","Disease association","https://platform-docs.opentargets.org/data-access/graphql-api","GraphQL v4","Ensembl target → associations"),
        ("RCSB PDB","Structure metadata","https://data.rcsb.org/","Data API v1","accession → structure metadata"),
        ("ChEMBL","Known chemistry and assays","https://www.ebi.ac.uk/chembl/api/data/docs","REST data API","target → assay measurements"),
        ("PubChem","Chemical identity","https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest","PUG REST","compound identifier → metadata"),
    ]
    rows=[]
    for name,role,doc,version,contract in source_data:
        rows.append(ToolSpec(name=name,role=role,documentation=doc,pinned_version=version,inputs=[contract.split(' → ')[0]],outputs=[contract.split(' → ')[1]],dependencies=['httpx==0.28.1'],external_data="Public target query / identifiers; disabled in local-only projects",cost_class="public API",implementation_status="implemented",validation_status="contract tests; see docs/retrieval.md for live verification",health_check="Requires network; cached responses remain usable").model_dump())
    local=[('RDKit','rdkit','2025.3.6','SMILES / SDF','canonical SDF + descriptors','https://www.rdkit.org/docs/'),('Meeko','meeko','0.7.1','prepared explicit-H SDF','ligand PDBQT','https://meeko.readthedocs.io/'),('AutoDock Vina','vina','1.2.7','receptor / ligand PDBQT + box','poses + Vina score (kcal/mol)','https://autodock-vina.readthedocs.io/'),('PoseBusters','posebusters','0.5.1','SDF pose + original ligand + receptor','physical validity flags','https://github.com/maabuu/posebusters'),('PLIP','plip','2.4.0','complex PDB','interaction XML','https://github.com/pharmai/plip'),('PDBFixer','pdbfixer','1.12','protein PDB','proposed repairs / prepared PDB','https://github.com/openmm/pdbfixer'),('PDB2PQR / PROPKA','pdb2pqr','3.7.1','supported protein PDB + explicit pH','PQR + protonation log','https://pdb2pqr.readthedocs.io/')]
    for name,module,version,inp,out,doc in local:
        available=importlib.util.find_spec(module) is not None
        if module=='vina': available=available or bool(shutil.which('vina')) or (Path(__file__).parent.parent/'.tools/vina').exists()
        if module=='plip': available=available or bool(shutil.which('plip'))
        rows.append(ToolSpec(name=name,role="Scientific computation",documentation=doc,pinned_version=version,inputs=[inp],outputs=[out],dependencies=[module],implementation_status="implemented supported adapter",validation_status="See docs/science.md; health does not establish scientific validity",health_check="available" if available else "missing dependency",code_terms="See pinned upstream license; PLIP GPL-2.0, Vina Apache-2.0, RDKit BSD-3-Clause").model_dump())
    rows.append(ToolSpec(name='OpenAI Responses',role='Read-only evidence advisor',documentation='https://developers.openai.com/api/docs/guides/function-calling',pinned_version='openai==2.26.0',inputs=['public project context + untrusted evidence'],outputs=['strict structured advisory proposal'],dependencies=['OPENAI_API_KEY'],external_data='Explicit opt-in; project context and cited evidence sent to OpenAI',cost_class='paid API',implementation_status='implemented',validation_status='mock-tested; live model unverified without credentials',health_check='configured' if os.getenv('OPENAI_API_KEY') else 'missing credential').model_dump())
    from .engines import engine_catalog
    for engine in engine_catalog():
        if engine['id'] == 'vina':
            continue  # Vina retains its established preparation/pose-check adapter above.
        row = ToolSpec(name=engine['name'], role=engine['role'], documentation=engine['docs'][0],
                       pinned_version='Runtime executable/model hashes recorded; local install required',
                       inputs=list(engine['input_contract']['properties']),
                       outputs=[f"{metric['name']} ({metric['unit']}; {metric['direction']})" for metric in engine['output_metrics']],
                       dependencies=[engine['id']], implementation_status=engine['implementation_status'],
                       validation_status=engine['validation_status'],
                       health_check='configured; validate installation before scientific use' if engine['available'] else '; '.join(engine['blockers'])).model_dump()
        row.update(engine_id=engine['id'], available=engine['available'], hardware=engine['hardware'])
        rows.append(row)
    for name in ('Boltz-2','AlphaFold prediction','ADMET-AI','OpenMM / GROMACS','OpenFE','Generative chemistry'):
        rows.append({'name':name,'implementation_status':'deferred','validation_status':'not implemented','health_check':'unavailable','role':'Optional later workflow'})
    return rows
