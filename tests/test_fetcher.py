import pytest
from unittest.mock import patch, MagicMock
from rdkit import Chem
from chemease import get_mol_from_cid, get_mol_from_name
from chemease.exceptions import FetchError, InvalidMoleculeError

@patch('chemease.fetcher.pcp.get_compounds')
def test_get_mol_from_name(mock_get_compounds):
    mock_compound = MagicMock()
    mock_compound.isomeric_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    mock_compound.cid = 2244
    mock_get_compounds.return_value = [mock_compound]

    mol = get_mol_from_name("aspirin")
    assert mol is not None
    assert mol.GetNumAtoms() > 0
    assert mol.HasProp("_Name")
    assert mol.GetProp("_Name") == "aspirin"
    mock_get_compounds.assert_called_once_with("aspirin", 'name')

@patch('chemease.fetcher.pcp.Compound.from_cid')
def test_get_mol_from_cid(mock_from_cid):
    mock_compound = MagicMock()
    mock_compound.isomeric_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    mock_compound.cid = 2244
    mock_compound.iupac_name = "2-acetoxybenzoic acid"
    mock_from_cid.return_value = mock_compound

    mol = get_mol_from_cid(2244)
    assert mol is not None
    assert mol.GetNumAtoms() > 0
    assert mol.GetProp("CID") == "2244"
    mock_from_cid.assert_called_once_with(2244)

@patch('chemease.fetcher.pcp.get_compounds')
def test_invalid_name(mock_get_compounds):
    mock_get_compounds.return_value = []
    with pytest.raises(FetchError):
        get_mol_from_name("thisisnotarealchemical123456789")
