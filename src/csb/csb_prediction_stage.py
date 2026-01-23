#!/usr/bin/python
from src.csb import csb_cluster
from src.db import database


def csb_finder(options) -> None:
    """
    Runs the CSB finder algorithm and updates the database with syntenic block clusters.

    Args:
        options (Options): Pipeline options

    Input Example:
        options.database_directory: str = "./results/my_db.sqlite"

    Output:
        Database is updated with new CSB clusters.
    """

    csb_cluster.csb_prediction(options)
    csb_gene_cluster_dict = csb_cluster.csb_jaccard(
        options, 0.0
    )  # 0.0 does merge csb but create the csb cluster dict

    database.index_database(options.database_directory)
    database.delete_keywords_from_csb(
        options.database_directory, options
    )  # remove keys with options.csb_name_prefix options.csb_name_suffix to avoid old keyword interference
    database.update_keywords(
        options.database_directory, csb_gene_cluster_dict
    )  # assigns the names of the keywords to the clusters
