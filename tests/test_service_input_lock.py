"""Partial comparisons freeze scientific inputs before a selected pilot exists."""
import pytest

from biocentric.service import Service


@pytest.mark.parametrize('record', ['comparison', 'selection_lock'])
@pytest.mark.parametrize('operation', ['upload', 'measurements'])
def test_comparison_or_lock_blocks_input_changes_during_revision(tmp_path, record, operation):
    service = Service(tmp_path)
    pid = service.create_project({'name': 'Frozen comparison', 'target': 'P00520',
                                  'organism': 'Mus musculus', 'workflow_mode': 'research'})
    service.store.research_update(pid, **{record: {'status': 'partial'},
                                 'input_revision': {'status': 'awaiting_inputs'}})
    with pytest.raises(ValueError, match='frozen|pilot|new research|new conversation'):
        if operation == 'upload':
            service.upload(pid, 'new.smi', b'CCO ethanol\n', 'library')
        else:
            service.import_measurements(pid, b'compound_id,endpoint,value,units\nc1,Ki,1,nM\n')
