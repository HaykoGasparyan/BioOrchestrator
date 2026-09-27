import copy
import json

import pytest

from biocentric.biocentric_model import (BioCentricInputs, BioCentricParameters, ENGINE_ID, MODELS,
                                       SUPPORTED_WEIGHTS_SHA256, parse_output, plan)


def payload():
    pair = {"record_index": 1, "compound_id": "ligand_A", "protein_id": "ABL1_chainA", "smiles": "CCO"}
    row = {**pair, "status": "completed", "score": 6.0, "per_model": {f"{name}@{epoch}": 6.0 for name, epoch in MODELS},
           "checkpoint_sha256": SUPPORTED_WEIGHTS_SHA256, "failure_reason": None, "ensemble_spread": 0.0}
    return {"schema_version": 1, "engine_id": ENGINE_ID, "expected_pairs": [pair], "predictions": [row]}


def test_native_scores_keep_identity_and_score_semantics(tmp_path):
    data = payload()
    path = tmp_path / "predictions.json"
    path.write_text(json.dumps(data))
    rows = parse_output(path, data["expected_pairs"])
    assert rows[0]["name"] == "ligand_A"
    assert rows[0]["metrics"][0]["unit"] == "pK"
    assert not rows[0]["metrics"][0]["experimental_measurement"]
    assert "pose_path" not in rows[0]


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "substitute", "nan", "boolean", "mean", "checkpoint", "model", "spread"])
def test_parser_rejects_incomplete_or_mislabelled_predictions(tmp_path, mutation):
    data = payload()
    expected = copy.deepcopy(data["expected_pairs"])
    row = data["predictions"][0]
    if mutation == "missing":
        data["predictions"] = []
    elif mutation == "duplicate":
        data["predictions"].append(copy.deepcopy(row))
    elif mutation == "substitute":
        row["protein_id"] = "other_protein"
    elif mutation == "nan":
        row["score"] = float("nan")
    elif mutation == "boolean":
        row["score"] = True
    elif mutation == "mean":
        row["score"] = 7.0
    elif mutation == "checkpoint":
        row["checkpoint_sha256"] = "0" * 64
    elif mutation == "model":
        row["per_model"]["unreviewed_model"] = row["per_model"].pop(next(iter(row["per_model"])))
    elif mutation == "spread":
        row["ensemble_spread"] = float("inf")
    path = tmp_path / "predictions.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        parse_output(path, expected)


def test_failed_pairs_are_retained_without_zero_or_fabricated_score(tmp_path):
    data = payload()
    data["predictions"][0].update(status="failed", score=None, per_model={}, failure_reason="Chemistry changed during reference PDB conversion")
    path = tmp_path / "predictions.json"
    path.write_text(json.dumps(data))
    row = parse_output(path, data["expected_pairs"])[0]
    assert row["metrics"] == []
    assert row["score"] is None
    assert row["failure_reason"]


def test_standalone_parser_bounds_native_output_before_json_decoding(tmp_path):
    path = tmp_path / 'predictions.json'
    path.write_bytes(b' ' * (16 * 1024 * 1024 + 1))
    with pytest.raises(ValueError, match='16 MiB'):
        parse_output(path, payload()['expected_pairs'])


def test_runtime_and_input_contract_are_bounded():
    with pytest.raises(ValueError):
        BioCentricParameters(device="automatically_pick_gpu")
    with pytest.raises(ValueError):
        BioCentricParameters(max_pairs=1001)
    with pytest.raises(ValueError):
        BioCentricInputs(receptor="r.pdb", ligands="l.sdf", protein_id="p", preparation_reviewed=False)


def _inputs(tmp_path, names=("one", "two"), smiles=("CCO", "CCN")):
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem
    receptor = tmp_path / "receptor.pdb"
    receptor.write_text("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 10.00           C  \nEND\n")
    ligands = tmp_path / "ligands.sdf"
    writer = Chem.SDWriter(str(ligands))
    for name, smi in zip(names, smiles):
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        mol.SetProp("_Name", name)
        writer.write(mol)
    writer.close()
    return {"receptor": str(receptor), "ligands": str(ligands), "protein_id": "test", "preparation_reviewed": True}


@pytest.mark.parametrize("names,smiles", [(("one", "one"), ("CCO", "CCN")), (("one", "two"), ("CCO", "CCO")), (("", "two"), ("CCO", "CCN"))])
def test_plan_rejects_duplicate_or_unnamed_pairs(tmp_path, monkeypatch, names, smiles):
    monkeypatch.delenv("BIOCENTRIC_MODEL_DIR", raising=False)
    with pytest.raises(ValueError, match="compound IDs|duplicate chemical"):
        plan(_inputs(tmp_path, names, smiles), {}, tmp_path / "out")


def test_absent_runtime_blocks_without_importing_torch(tmp_path, monkeypatch):
    monkeypatch.delenv("BIOCENTRIC_MODEL_DIR", raising=False)
    monkeypatch.delenv("BIOCENTRIC_MODEL_PYTHON", raising=False)
    result = plan(_inputs(tmp_path), {}, tmp_path / "out")
    assert result["status"] == "blocked"
    assert result["expected_pairs"][0]["compound_id"] == "one"
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize('directory_state', ['missing', 'empty', 'interrupted'])
def test_resume_plan_revalidates_exact_request_without_changing_evidence(tmp_path, monkeypatch, directory_state):
    from biocentric import engines
    monkeypatch.delenv('BIOCENTRIC_MODEL_DIR', raising=False)
    monkeypatch.delenv('BIOCENTRIC_MODEL_PYTHON', raising=False)
    inputs = _inputs(tmp_path)
    output = tmp_path / 'out'
    reviewed = engines.plan_engine_run(ENGINE_ID, inputs, {}, output)
    if directory_state != 'missing':
        output.mkdir()
    if directory_state == 'interrupted':
        (output / 'request.json').write_text(reviewed['generated_files']['request.json'])
        (output / 'engine.log').write_text('Preserved interrupted native log')
    before = {path.name: path.read_bytes() for path in output.iterdir()} if output.exists() else {}
    resumed = engines.plan_engine_run(ENGINE_ID, inputs, {}, output, resume_plan=reviewed)
    assert resumed == reviewed
    assert ({path.name: path.read_bytes() for path in output.iterdir()} if output.exists() else {}) == before


@pytest.mark.parametrize('change', ['request', 'missing_request', 'escaped_request', 'oversize_request', 'parameters',
                                  'receptor', 'command', 'adapter', 'runtime'])
def test_resume_rejects_request_input_protocol_or_runtime_substitution(tmp_path, monkeypatch, change):
    from pathlib import Path
    from biocentric import biocentric_model
    monkeypatch.delenv('BIOCENTRIC_MODEL_DIR', raising=False)
    monkeypatch.delenv('BIOCENTRIC_MODEL_PYTHON', raising=False)
    inputs = _inputs(tmp_path)
    output = tmp_path / 'out'
    reviewed = plan(inputs, {}, output)
    output.mkdir()
    request_path = output / 'request.json'
    request_path.write_text(reviewed['generated_files']['request.json'])
    (output / 'engine.log').write_text('Original execution evidence')
    parameters = {}
    if change == 'request':
        request_path.write_text('{}')
    elif change == 'missing_request':
        request_path.unlink()
    elif change == 'escaped_request':
        outside = tmp_path / 'outside.json'
        outside.write_bytes(request_path.read_bytes())
        request_path.unlink()
        request_path.symlink_to(outside)
    elif change == 'oversize_request':
        request_path.write_bytes(b' ' * (16 * 1024 * 1024 + 1))
    elif change == 'parameters':
        parameters = {'cpu': 1}
    elif change == 'receptor':
        receptor = Path(inputs['receptor'])
        receptor.write_text(receptor.read_text() + 'REMARK Changed evidence\n')
    elif change == 'command':
        reviewed['commands'][0][0] = 'unreviewed-interpreter'
    elif change == 'adapter':
        reviewed['provenance']['adapter_sha256'] = 'f' * 64
    else:
        installation = biocentric_model.installation()
        monkeypatch.setattr(biocentric_model, 'installation', lambda: {**installation, 'code_sha256': 'other-code'})
    with pytest.raises((ValueError, FileNotFoundError)):
        plan(inputs, parameters, output, resume_plan=reviewed)
    assert (output / 'engine.log').read_text() == 'Original execution evidence'


def test_other_native_engines_do_not_accept_biocentric_resume_plans(tmp_path):
    from biocentric import engines
    with pytest.raises(ValueError, match='only supported for BioCentric'):
        engines.plan_engine_run('vina', {}, {}, tmp_path / 'out', resume_plan={})


def test_receptor_heterogens_require_explicit_preparation(tmp_path):
    inputs = _inputs(tmp_path)
    from pathlib import Path
    path = Path(inputs["receptor"])
    path.write_text(path.read_text().replace("ATOM  ", "HETATM"))
    with pytest.raises(ValueError, match="ATOM records"):
        plan(inputs, {}, tmp_path / "out")


@pytest.mark.parametrize("failure_count,status", [(0, "completed"), (1, "partial"), (2, "failed")])
def test_outer_engine_execution_accounts_for_every_pair(tmp_path, monkeypatch, failure_count, status):
    """Fixture subprocess contract only; real inference has a separate verifier."""
    from biocentric import biocentric_model, engines, science
    original = biocentric_model.installation
    monkeypatch.delenv("BIOCENTRIC_MODEL_DIR", raising=False)
    monkeypatch.delenv("BIOCENTRIC_MODEL_PYTHON", raising=False)
    monkeypatch.setattr(biocentric_model, "installation", lambda: {**original(), "blockers": [], "weights_sha256": SUPPORTED_WEIGHTS_SHA256})
    def runner(argv, cwd, log_path, timeout, cancelled):
        from pathlib import Path
        request = json.loads(Path(argv[-1]).read_text())
        rows = []
        for index, pair in enumerate(request["pairs"]):
            row = {**payload()["predictions"][0], **pair, 'code_sha256': request['code_sha256'],
                   'device': request['parameters']['device'], 'models': request['models']}
            if index < failure_count:
                row.update(status="failed", score=None, per_model={}, failure_reason="Synthetic unsupported preparation fixture")
            rows.append(row)
        (cwd / "predictions.json").write_text(json.dumps({"schema_version": 1, "engine_id": ENGINE_ID,
            "expected_pairs": request["pairs"], "predictions": rows}))
        log_path.write_text("Synthetic subprocess fixture; no model executed\n")
        return {"returncode": 0}
    monkeypatch.setattr(science, "run_process", runner)
    result = engines.run_engine(ENGINE_ID, _inputs(tmp_path), {}, tmp_path / "out")
    assert result["status"] == status
    assert result["accounting"] == {"requested": 2, "completed": 2 - failure_count, "failed": failure_count}
    assert len(result["results"]) == 2
    assert sum(not row["metrics"] for row in result["results"]) == failure_count


def test_checkpoint_reuse_verifies_request_content_and_exact_pair_identity(tmp_path):
    from biocentric.biocentric_predictor import read_checkpoint, row_sha256
    row = payload()["predictions"][0]
    expected = payload()["expected_pairs"][0]
    path = tmp_path / "000001.json"
    checkpoint = {"request_sha256": "a" * 64, "row_sha256": row_sha256(row), "row": row}
    path.write_text(json.dumps(checkpoint))
    assert read_checkpoint(path, "a" * 64, expected) == row
    with pytest.raises(ValueError, match="hash verification"):
        read_checkpoint(path, "b" * 64, expected)
    original = path.read_bytes()
    with pytest.raises(ValueError, match="identity"):
        read_checkpoint(path, "a" * 64, {**expected, "protein_id": "other"})
    assert path.read_bytes() == original
    checkpoint["row"]["score"] = 42.0
    path.write_text(json.dumps(checkpoint))
    changed = path.read_bytes()
    with pytest.raises(ValueError, match="hash verification"):
        read_checkpoint(path, "a" * 64, expected)
    assert path.read_bytes() == changed  # Original invalid evidence is never overwritten.


def test_oom_recovery_releases_unreachable_objects_before_gpu_cache_cleanup():
    """A synthetic allocation cycle tests cleanup order; this is not a GPU test."""
    import gc
    from types import SimpleNamespace
    import weakref
    from biocentric.biocentric_predictor import recover_cuda_after_failure
    class Allocation:
        pass
    allocation = Allocation()
    allocation.cycle = allocation
    reference = weakref.ref(allocation)
    calls = []
    def empty_cache():
        assert reference() is None
        calls.append("released")
    enabled = gc.isenabled()
    gc.disable()
    try:
        del allocation
        assert reference() is not None
        recover_cuda_after_failure(SimpleNamespace(cuda=SimpleNamespace(empty_cache=empty_cache)))
        assert calls == ["released"]
    finally:
        if enabled:
            gc.enable()


def test_explicit_pose_mode_allows_repeated_chemistry_but_not_repeated_ids(tmp_path, monkeypatch):
    monkeypatch.delenv("BIOCENTRIC_MODEL_DIR", raising=False)
    inputs = _inputs(tmp_path, smiles=("CCO", "CCO"))
    result = plan(inputs, {"input_mode": "pose_records"}, tmp_path / "out")
    assert len(result["expected_pairs"]) == 2
    assert result["parameters"]["input_mode"] == "pose_records"
    with pytest.raises(ValueError, match="compound IDs"):
        plan(_inputs(tmp_path, names=("one", "one")), {"input_mode": "pose_records"}, tmp_path / "other")


def test_pdb_missing_double_bond_labels_are_verified_from_geometry():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from biocentric.biocentric_predictor import canonical, validate_pdb_identity
    original = Chem.AddHs(Chem.MolFromSmiles('C/C=C/C'))
    assert AllChem.EmbedMolecule(original, randomSeed=42) == 0
    roundtrip = Chem.MolFromPDBBlock(Chem.MolToPDBBlock(original), removeHs=False)
    Chem.RemoveStereochemistry(roundtrip)  # PDB carries coordinates but no E/Z field.
    before = roundtrip.GetConformer().GetPositions().copy()
    assert validate_pdb_identity(roundtrip, canonical(original)) == 'verified_3d_stereochemistry'
    assert (roundtrip.GetConformer().GetPositions() == before).all()
    with pytest.raises(ValueError, match='stereochemistry'):
        validate_pdb_identity(roundtrip, 'C/C=C\\C')
    changed = Chem.Mol(roundtrip)
    changed.GetAtomWithIdx(0).SetFormalCharge(1)
    with pytest.raises(ValueError):
        validate_pdb_identity(changed, canonical(original))
