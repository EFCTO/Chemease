from rdkit import Chem
from rdkit.Chem import AllChem


def run_virtual_reaction(
    rxn_smarts: str, reactants: list[Chem.Mol]
) -> list[list[Chem.Mol]]:
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    if rxn is None:
        raise ValueError(f"Invalid reaction SMARTS: {rxn_smarts}")

    products = rxn.RunReactants(tuple(reactants))
    return products


# Common reaction definitions
AMIDE_COUPLING_SMARTS = (
    "[CX3](=[OX1])[F,Cl,Br,I,OH1,O-].[NX3;H2,H1:1]>>[CX3](=[OX1])[N:1]"
)
CLICK_CHEMISTRY_SMARTS = "[C:1]#[C:2].[N:3]=[N+:4]=[N-:5]>>[c:1]1[c:2][n:3][n:4][n:5]1"


def synthesize_amides(acids: list[Chem.Mol], amines: list[Chem.Mol]) -> list[Chem.Mol]:
    rxn = AllChem.ReactionFromSmarts(AMIDE_COUPLING_SMARTS)
    results = []
    for acid in acids:
        if acid is None:
            continue
        for amine in amines:
            if amine is None:
                continue
            products = rxn.RunReactants((acid, amine))
            for prod_set in products:
                # Amide coupling produces one main product usuall
                mol = prod_set[0]
                Chem.SanitizeMol(mol)
                results.append(mol)
    return results

__v