import click
import pandas as pd
from rdkit import Chem

from chemease import get_mol_from_name
from chemease.analysis import cluster_molecules
from chemease.descriptors import calculate_advanced_descriptors_df
from chemease.io import read_sdf, write_sdf
from chemease.logger import logger


@click.group()
def cli():
    pass


@cli.command()
@click.option("--name", required=True, help="Name of the molecule to fetch.")
@click.option("--output", "-o", required=True, help="Output SDF file path.")
def fetch(name: str, output: str):
    try:
        mol = get_mol_from_name(name)
        write_sdf([mol], output)
        logger.info(f"Successfully fetched {name} and saved to {output}")
    except Exception as e:
        logger.error(f"Error fetching {name}: {e}")


@cli.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option("--output", "-o", required=True, help="Output CSV file path.")
def descriptors(input_file: str, output: str):
    try:
        mols = read_sdf(input_file)
        logger.info(f"Read {len(mols)} molecules from {input_file}")
        df = calculate_advanced_descriptors_df(mols)
        df.to_csv(output, index=False)
        logger.info(f"Calculated descriptors and saved to {output}")
    except Exception as e:
        logger.error(f"Error processing {input_file}: {e}")


@cli.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "--cutoff", default=0.5, help="Clustering distance cutoff (default: 0.5)."
)
def cluster(input_file: str, cutoff: float):
    try:
        mols = read_sdf(input_file)
        logger.info(f"Clustering {len(mols)} molecules with cutoff {cutoff}...")
        clusters = cluster_molecules(mols, cutoff=cutoff)
        logger.info(f"Found {len(clusters)} clusters.")
        for i, c in enumerate(clusters[:5]):
            logger.info(f"Cluster {i+1} size: {len(c)}")
        if len(clusters) > 5:
            logger.info("...")
    except Exception as e:
        logger.error(f"Error clustering {input_file}: {e}")


if __name__ == "__main__":
    cli()
