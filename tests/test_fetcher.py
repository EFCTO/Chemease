import pytest

from chemease import get_mol_from_cid, get_mol_from_name


def test_get_mol_from_name():
    mol = get_mol_from_name("aspirin")
    assert mol is not None
    assert mol.GetNumAtoms() > 0
    assert mol.HasProp("_Name")
    assert mol.GetProp("_Name") == "aspirin"


def test_get_mol_from_cid():
    mol = get_mol_from_cid(2244)  # Aspirin CID
    assert mol is not None
    assert mol.GetNumAtoms() > 0
    assert mol.GetProp("CID") == "2244"


def test_invalid_name():
    with pytest.raises(Exception):
        get_mol_from_name("thisisnotarealchemical123456789")
