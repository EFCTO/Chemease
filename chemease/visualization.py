from rdkit import Chem
from rdkit.Chem import Draw


def draw_molecule(mol: Chem.Mol, size: tuple[int, int] = (300, 300)):
    # Draw a single molecule Returns a PIL Image
    if mol is None:
        return None
    return Draw.MolToImage(mol, size=size)


def draw_grid(
    mols: list[Chem.Mol],
    molsPerRow: int = 3,
    subImgSize: tuple[int, int] = (200, 200),
    legends: list[str] = None,
):
    # Draw a grid of molecules Returns a PIL Image
    valid_mols = []
    valid_legends = []
    for i, mol in enumerate(mols):
        if mol is not None:
            valid_mols.append(mol)
            if legends and i < len(legends):
                valid_legends.append(legends[i])
            else:
                valid_legends.append(
                    mol.GetProp("_Name") if mol.HasProp("_Name") else ""
                )

    if not valid_mols:
        return None

    return Draw.MolsToGridImage(
        valid_mols, molsPerRow=molsPerRow, subImgSize=subImgSize, legends=valid_legends
    )

def i_like_caffeine():
    return "Caffeine is great!"
