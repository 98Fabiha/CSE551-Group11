import pandas as pd
import itertools
from itertools import combinations
from math import factorial

def find_unique_genes(dataframe):
    """
    Identify unique genes for each cancer type.

    Parameters:
        dataframe (pd.DataFrame): The dataset containing cancer types and gene data.

    Returns:
        dict: A dictionary mapping each cancer type to its unique genes.
    """
    unique_genes = {}
    all_cancers = dataframe['Cancer'].unique()

    for cancer in all_cancers:
        # Extract gene data for the current cancer type
        cancer_row = dataframe.loc[dataframe['Cancer'] == cancer].iloc[0, 1:]
        cancer_genes = set(cancer_row[cancer_row == 1].index)

        # Extract gene data for all other cancers
        other_rows = dataframe.loc[dataframe['Cancer'] != cancer].iloc[:, 1:]
        other_genes = set(other_rows.columns[(other_rows.sum(axis=0) > 0)])

        # Find genes unique to the current cancer
        unique_genes[cancer] = cancer_genes - other_genes

    return unique_genes


def find_minimal_unique_subsets(unique_genes):
    """
    Calculate minimal unique subsets for each cancer type.

    Parameters:
        unique_genes (dict): A dictionary of unique genes per cancer type.

    Returns:
        dict: A dictionary mapping each cancer type to its minimal unique subsets.
    """
    minimal_subsets = {}

    for cancer, genes in unique_genes.items():
        if not genes:
            minimal_subsets[cancer] = []  # No unique genes
            continue

        gene_list = list(genes)
        minimal_set = []

        # Check all combinations of genes
        for size in range(1, len(gene_list) + 1):
            for subset in combinations(gene_list, size):
                subset_set = set(subset)

                # Check if this subset uniquely identifies the cancer type
                if all(
                    subset_set.isdisjoint(other)  # The subset must not be shared with other cancers
                    for other in minimal_set
                ):
                    minimal_set.append(subset_set)

        minimal_subsets[cancer] = minimal_set  # Store all valid minimal subsets

    return minimal_subsets


def calculate_fk_msu_fm_suda(unique_genes, minimal_subsets):
    """
    Calculate FK, MSU, fM, and SUDA values for each cancer type.

    Parameters:
        unique_genes (dict): A dictionary of unique genes per cancer type.
        minimal_subsets (dict): A dictionary of minimal unique subsets per cancer type.

    Returns:
        tuple: Four dictionaries containing FK, MSU, fM, and SUDA values for each cancer type.
    """
    fk_values = {}
    msu_values = {}
    fm_values = {}
    suda_values = {}

    for cancer, subsets in minimal_subsets.items():
        # FK: Frequency of unique genes
        fk_values[cancer] = len(unique_genes[cancer])

        # fM: Frequency of minimal unique subsets (count of minimal subsets)
        fm_values[cancer] = len(subsets)

        # MSU: The value is 1 because there are multiple minimal subsets
        msu_values[cancer] = 1 if minimal_subsets[cancer] else 0  # MSU will be 1 regardless of number of subsets

        # SUDA: Sum of factorials for the number of attributes in each subset
        suda_score = sum(factorial(len(subset)) for subset in subsets)
        suda_values[cancer] = suda_score

    return fk_values, msu_values, fm_values, suda_values

def calculate_all_subsets(unique_genes, sample_size=1000):
    """
    Calculate all possible subsets for each cancer type and sample a number of them.

    Parameters:
        unique_genes (dict): A dictionary of unique genes per cancer type.
        sample_size (int): Number of subsets to sample for each cancer type.

    Returns:
        dict: A dictionary mapping each cancer type to its sampled subsets.
    """
    all_subsets = {}
    for cancer, genes in unique_genes.items():
        if not genes:
            all_subsets[cancer] = []
            continue

        gene_list = list(genes)
        all_combinations = []

        # Generate combinations one size at a time
        for size in range(1, len(gene_list) + 1):
            combinations = list(itertools.combinations(gene_list, size))
            all_combinations.extend(combinations)

        # Limit the number of combinations to the sample_size
        if len(all_combinations) > sample_size:
            all_combinations = all_combinations[:sample_size]

        all_subsets[cancer] = all_combinations

    return all_subsets

