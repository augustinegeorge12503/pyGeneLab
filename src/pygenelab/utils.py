# utils.py 

"""
all utility functions
"""

# imports
from pathlib import Path
from itertools import chain, repeat

import pandas as pd
from scipy import stats


# convert_gmt_to_decoupler_format
def convert_gmt_to_decoupler_format(
    pth: Path,
    include_pathways=None,
    gene_origin="mice"
) -> pd.DataFrame:
    """
    convert .gmt file paths to decoupler input format
    """

    # convert_gmt_to_decoupler_format
    # api:
    # convert_gmt_to_decoupler_format(
    #     pth=gmt_path,
    #     include_pathways=["PATHWAY_1", "PATHWAY_2"],
    #     gene_origin="mice"
    # )

    # check gene origin
    if gene_origin not in ["mice", "human"]:
        raise ValueError("gene_origin must be either 'mice' or 'human'")

    # make pathway filter set
    if include_pathways is not None:
        include_pathways = set(include_pathways)

    # dictionary to store selected pathways
    pathways = {}

    # open .gmt path and get pathway: genes
    with Path(pth).open("r") as f:
        for line in f:
            name, _, *genes = line.strip().split("\t")

            # skip pathways not in selected list
            if include_pathways is not None and name not in include_pathways:
                continue

            # format gene names
            if gene_origin == "mice":
                genes = [gene.capitalize() for gene in genes]

            elif gene_origin == "human":
                genes = [gene.upper() for gene in genes]

            pathways[name] = genes

    # decoupler accepts "source" for pathway and "target" for genes
    return pd.DataFrame.from_records(
        chain.from_iterable(zip(repeat(k), v) for k, v in pathways.items()),
        columns=["source", "target"],
    )


# calculate_pairwise_significance
def calculate_pairwise_significance(data, groups, x_var, y_var):
    """
    calculate pairwise mann-whitney significance between groups
    """

    results = {}

    # compare every pair of groups
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):

            # get values for each group
            group1 = data[data[x_var] == groups[i]][y_var].dropna()
            group2 = data[data[x_var] == groups[j]][y_var].dropna()

            # skip if one group is empty
            if len(group1) == 0 or len(group2) == 0:
                results[(i, j)] = {
                    "p-value": None,
                    "significance": "ns"
                }
                continue

            # run mann-whitney u test
            statistic, pvalue = stats.mannwhitneyu(
                group1,
                group2,
                alternative="two-sided"
            )

            # assign significance stars
            if pvalue < 0.001:
                sig = "***"
            elif pvalue < 0.01:
                sig = "**"
            elif pvalue < 0.05:
                sig = "*"
            else:
                sig = "ns"

            # store result using group positions
            results[(i, j)] = {
                "p-value": pvalue,
                "significance": sig
            }

    # return pairwise results
    return results
