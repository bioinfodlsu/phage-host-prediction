"""
=========================================================================
This script contains the utility functions for exploratory data analysis.
@author    Mark Edward M. Gonzales
=========================================================================
"""

import os


class EDAUtil(object):
    def __init__(self):
        """
        Constructor
        """
        pass

    def get_rbps(self, directory):
        """
        Get the protein IDs in the RBP embeddings files stored in the given directory

        Parameters:
        - directory: File path of the directory where the RBP embeddings files are saved

        Returns:
        - List of protein IDs in the RBP embeddings files stored in the directory
        """
        rbps = []
        for file in os.listdir(f"{directory}"):
            with open(f"{directory}/{file}") as embeddings:
                for embedding in embeddings:
                    rbp_id = embedding.split(",")[0]
                    if rbp_id != "ID":
                        rbps.append(rbp_id)

        return rbps
