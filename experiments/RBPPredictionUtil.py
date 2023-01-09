# ================================================================================================================
# This script contains the utility functions for the computational prediction of receptor-binding proteins (RBPs).
# ================================================================================================================

import os
import shutil
from collections import defaultdict

from SequenceParsingUtil import SequenceParsingUtil

import pandas as pd
import numpy as np

from xgboost import XGBClassifier

class RBPPredictionUtil(object):
    def __init__(self):
        pass
    
    def predict_rbps(self, xgb_rbp_prediction, hypothetical_genbank_dir, complete_genbank_dir):
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
                    # 1024 is the size of the dense vector representation yielded by ProtTransBert
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
        util = SequenceParsingUtil()
        
        with open(f'{prottransbert_complete_dir}/{phage}-rbp-embeddings.csv', 'r') as embeddings:
            rbp_prottransbert, _ = util.get_proteins_embeddings(embeddings)
            
        for protein in rbp_prottransbert:
            entry = f'>{protein}\n{self.get_sequence(protein, phage, True, *fasta_folders)}\n'
            
            with open(f'{dest_dir}/{phage}-rbp.fasta', 'a') as fasta:
                fasta.write(entry)