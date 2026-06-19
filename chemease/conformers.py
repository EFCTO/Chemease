from rdkit import Chem
from rdkit.Chem import AllChem


def generate_conformers(
    mol: Chem.Mol, num_confs: int = 10, optimize: bool = True
) -> Chem.Mol:
    """
    generate 3D conformers for a molecule using ETKDG.
    optionally optimizes the generated conformers using MMFF94.
    returns a new molecule object with explicit hydrogens and 3D conformers.
    """
    if mol is None:
        return None

    # Add explicit hydrogens for 3D generation
    mol_3d = Chem.AddHs(mol)

    # Generate conformers using ETKDG
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMultipleConfs(mol_3d, numConfs=num_confs, params=params)

    # Optimize conformers
    if optimize:
        AllChem.MMFFOptimizeMoleculeConfs(mol_3d, maxIters=500)

    return mol_3d


def get_conformer_energies(mol_3d: Chem.Mol) -> list[float]:
    if mol_3d is None or mol_3d.GetNumConformers() == 0:
        return []

    energies = []
    mp = AllChem.MMFFGetMoleculeProperties(mol_3d, mmffVariant="MMFF94")
    if mp is None:
        return []  # Cannot get properties, perhaps missing parameters for some atoms

    for conf in mol_3d.GetConformers():
        ff = AllChem.MMFFGetMoleculeForceField(mol_3d, mp, confId=conf.GetId())
        if ff:
            energies.append(ff.CalcEnergy())
        else:
            energies.append(None)
    return energies
