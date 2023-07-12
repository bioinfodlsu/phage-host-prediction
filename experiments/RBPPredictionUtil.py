"""
================================================================================================================
This script contains the utility functions for the computational prediction of receptor-binding proteins (RBPs).
@author    Mark Edward M. Gonzales
================================================================================================================
"""

import os

from SequenceParsingUtil import SequenceParsingUtil

import pandas as pd
import numpy as np

from xgboost import XGBClassifier

class RBPPredictionUtil(object):
    def __init__(self):
        """
        Constructor
        """
        pass
    
    def predict_rbps(self, xgb_rbp_prediction, hypothetical_genbank_dir, complete_genbank_dir):
        """
        Write the embeddings of the predicted RBPs to the pertinent CSV files, with each CSV file corresponding to a phage
        
        Parameters:
        - xgb_rbp_prediction: File path of the XGBoost model for the computational prediction of RBPs
                              (taken from https://github.com/dimiboeckaerts/PhageRBPdetection/blob/main/data/RBPdetect_xgb_model.json)
        - hypothetical_genbank_dir: File path of the directory containing the ProtBert embeddings of the annotated
                                    hypothetical proteins
        - complete_genbank_dir: File path of the directory containing the ProtBert embeddings of the annotated RBPs,
                                along with those that will be computationally predicted as RBPs
        """
        for file in os.listdir(hypothetical_genbank_dir):
            embeddings_df = pd.read_csv(f'{hypothetical_genbank_dir}/{file}')
            embeddings = np.asarray(embeddings_df.iloc[:, 1:])
            protein_ids = np.asarray(embeddings_df.iloc[:, :1])

            xgb_saved = XGBClassifier()
            xgb_saved.load_model(xgb_rbp_prediction)

            score_xgb = xgb_saved.predict_proba(embeddings)[:, 1]
            preds_xgb = (score_xgb > 0.5) * 1

            # Write to file only if there is at least one predicted RBP
            if 1 in preds_xgb:        
                filename = file[:-len('-hypothetical-embeddings.csv')] + '-rbp-embeddings.csv'    
                filename = f'{complete_genbank_dir}/{filename}'

                print(file)

                # Add header if file is newly created
                newly_created = False
                if not os.path.exists(filename):
                    header = 'ID'
                    # 1024 is the size of the dense vector representation yielded by ProtBert
                    for i in range(1024):
                        header += ',' + str(i)
                    header += '\n'

                    newly_created = True

                with open(f'{filename}', 'a') as file:
                    if newly_created:
                        file.write(header)
                    for idx, pred in enumerate(preds_xgb):
                        if pred == 1:
                            embedding = ','.join(map(str, embeddings[idx]))
                            entry = protein_ids[idx][0] + ',' + embedding + '\n'
                            file.write(entry)
                            
    def check_differences(self, prottransbert_complete_dir, plm_complete_dir):
        """
        Prints the difference between the protein IDs in the ProtBert embeddings folder and in the embeddings folder
        of a given protein language model. There should be a one-to-one correspondence between the protein IDs in these
        two directories
        
        Parameters:
        - prottransbert_complete_dir: File path of the directory containing the ProtBert embeddings of the annotated RBPs,
                                      along with those that will be computationally predicted as RBPs
        - plm_complete_dir: File path of the directory containing the protein embeddings of the annotated RBPs,
                            along with those that will be computationally predicted as RBPs
        
        Returns:
        - Set of phages where there are differences in the protein IDs in the ProtBert embeddings folder and in the embeddings
          folder of the protein language model
        """
        erroneous = set()
        util = SequenceParsingUtil()

        for file in os.listdir(prottransbert_complete_dir):
            phage = file[:-len('-rbp-embeddings.csv')]

            with open(f'{prottransbert_complete_dir}/{file}', 'r') as embeddings:
                rbp_prottransbert, _ = util.get_proteins_embeddings(embeddings)

            try:
                with open(f'{plm_complete_dir}/{file}', 'r') as embeddings:
                    rbp_plm, _ = util.get_proteins_embeddings(embeddings)

                if rbp_prottransbert - rbp_plm or rbp_plm - rbp_prottransbert:
                    print(phage, '|', rbp_prottransbert - rbp_plm, '|', rbp_plm - rbp_prottransbert, '|')
                    erroneous.add(phage)

            except FileNotFoundError:
                print(phage, '| Missing')
                erroneous.add(phage)
                continue
                
        extra_files = set(os.listdir(plm_complete_dir)) - set(os.listdir(prottransbert_complete_dir))
        for file in extra_files:
            phage = file[:-len('-rbp-embeddings.csv')]
            print(phage, '| Extra')
            erroneous.add(phage)
                
        return erroneous
    
    def get_sequence(self, protein_id, phage_accession, protein, *fasta_folders):
        """
        Retrieves the sequence associated with a given protein ID
        
        Parameters:
        - protein_id: Protein ID
        - phage_accession: Accession ID of the phage to which the protein belongs
        - fasta_folders: Folders containing the FASTA files of the proteomes
        
        Returns:
        - Sequence associated with the protein ID
        """
        util = SequenceParsingUtil()
        
        extension = ''
        if protein:
            extension = 'fasta'
        else:
            extension = 'ffn'
        
        for folder in fasta_folders:
            filename = f'{folder}/{phage_accession}-rbp.{extension}'
            
            if os.path.exists(filename):
                sequence = util.get_seq_in_fasta(protein_id, filename)
                if sequence:
                    return sequence
                
            filename = f'{folder}/{phage_accession}-hypothetical.{extension}'
            
            if os.path.exists(filename):
                sequence = util.get_seq_in_fasta(protein_id, filename)
                if sequence:
                    return sequence
                
        return None
                
    def generate_fasta_addl_phages(self, phage, prottransbert_complete_dir, dest_dir, *fasta_folders):
        """
        Generate FASTA files containing the RBP sequences (including computationally predicted RBPs) of the given phage 
        
        - phage: Accession number of the phage
        - prottransbert_complete_dir: File path of the directory containing the ProtBert embeddings of the annotated RBPs,
                                      along with those that will be computationally predicted as RBPs
        - dest_dir: File path of the directory to which the FASTA files will be saved
        """
        util = SequenceParsingUtil()
        
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)
        
        with open(f'{prottransbert_complete_dir}/{phage}-rbp-embeddings.csv', 'r') as embeddings:
            rbp_prottransbert, _ = util.get_proteins_embeddings(embeddings)
            
        for protein in rbp_prottransbert:
            sequence = self.get_sequence(protein, phage, True, *fasta_folders)
            assert sequence != None
            
            entry = f'>{protein}\n{sequence}\n'
            
            with open(f'{dest_dir}/{phage}-rbp.fasta', 'a') as fasta:
                fasta.write(entry)