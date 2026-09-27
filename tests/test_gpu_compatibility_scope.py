"""Do not extrapolate the verified legacy Uni-Dock binary to untested inputs."""
import pytest

from biocentric.engines import Box, _review_unidock_110


@pytest.fixture
def inputs(tmp_path):
    receptor, ligand = tmp_path / 'receptor.pdbqt', tmp_path / 'ligand.pdbqt'
    receptor.write_text('ATOM      1 carbon C\n')
    ligand.write_text('REMARK SMILES CCO\nATOM      1 carbon C\nTORSDOF 0\n')
    return receptor, ligand


def test_legacy_unidock_common_chemistry_and_small_box_allowed(inputs):
    receptor, ligand = inputs
    _review_unidock_110(receptor, [ligand], Box(center=(0, 0, 0), size=(20, 20, 20)))


@pytest.mark.parametrize('side', [30, 35])
def test_legacy_unidock_large_box_rejected(inputs, side):
    receptor, ligand = inputs
    with pytest.raises(ValueError, match='below 30'):
        _review_unidock_110(receptor, [ligand], Box(center=(0, 0, 0), size=(side, 20, 20)))


@pytest.mark.parametrize('record', [
    'REMARK SMILES C1CCCCCCC1\nATOM      1 carbon C\n',
    'REMARK SMILES CCO\nATOM      1 glue CG0\n',
    'REMARK SMILES B\nATOM      1 boron B\n',
    'ATOM      1 carbon C\n',
])
def test_legacy_unidock_unsupported_or_unidentifiable_chemistry_rejected(inputs, record):
    receptor, ligand = inputs
    ligand.write_text(record)
    with pytest.raises(ValueError, match='unsupported|macrocycle|metadata'):
        _review_unidock_110(receptor, [ligand], Box(center=(0, 0, 0), size=(20, 20, 20)))
