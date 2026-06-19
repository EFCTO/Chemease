__version__ = "0.2.0"
__version__ = "0.1.0"

from .analysis import cluster_molecules, get_murcko_scaffold, perform_pca
from .conformers import generate_conformers, get_conformer_energies
from .databases import (get_bioactivities, get_mol_from_chembl,
                        get_target_by_name)
from .descriptors import calculate_advanced_descriptors_df, get_maccs_keys
from .fetcher import get_mol_from_cid, get_mol_from_name
from .filters import filter_by_substructure, passes_lipinski
from .io import read_sdf, read_smiles, write_sdf
from .properties import calculate_properties, calculate_properties_df
from .reactions import (AMIDE_COUPLING_SMARTS, CLICK_CHEMISTRY_SMARTS,
                        run_virtual_reaction, synthesize_amides)
from .similarity import calculate_similarity, find_similar
from .visualization import draw_grid, draw_molecule

__all__ = [
    "get_mol_from_name",
    "get_mol_from_cid",
    "calculate_properties",
    "calculate_properties_df",
    "passes_lipinski",
    "filter_by_substructure",
    "calculate_similarity",
    "find_similar",
    "draw_molecule",
    "draw_grid",
    "read_sdf",
    "write_sdf",
    "read_smiles",
    "get_target_by_name",
    "get_bioactivities",
    "get_mol_from_chembl",
    "get_maccs_keys",
    "calculate_advanced_descriptors_df",
    "generate_conformers",
    "get_conformer_energies",
    "get_murcko_scaffold",
    "cluster_molecules",
    "perform_pca",
    "run_virtual_reaction",
    "synthesize_amides",
    "AMIDE_COUPLING_SMARTS",
    "CLICK_CHEMISTRY_SMARTS",
]
