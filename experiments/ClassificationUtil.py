# ==============================================================
# This script contains the utility functions for classification.
# ==============================================================

import os
import copy

from RBPPredictionUtil import RBPPredictionUtil

import pandas as pd
import numpy as np

from ete3 import NCBITaxa
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, precision_recall_fscore_support
from Bio import SeqIO

class ClassificationUtil(object):
    def __init__(self, complete_embeddings_dir = None, RANDOM_NUM = 42):
        self.complete_embeddings_dir = complete_embeddings_dir
        self.RANDOM_NUM = RANDOM_NUM
        self.inphared_gb = None
    
    # =============
    # Preprocessing
    # =============
    def set_inphared_gb(self, inphared_gb):
        self.inphared_gb = inphared_gb
    
    def get_phages(self):
        phages = set()
        for file in os.listdir(f'{self.complete_embeddings_dir}'):
            phages.add(file[:-len('-rbp-embeddings.csv')])
            
        return phages
    
    def get_rbps(self):
        rbps = []
        for file in os.listdir(f'{self.complete_embeddings_dir}'):
            phage_name = file[:-len('-rbp-embeddings.csv')]
            with open(f'{self.complete_embeddings_dir}/{file}', 'r') as embeddings:
                for embedding in embeddings:
                    rbp_id = embedding.split(',')[0]
                    if rbp_id != 'ID':
                        rbps.append([rbp_id, phage_name])
                        
        return rbps
    
    def get_host_taxonomy(self, rbp_embeddings, host_column = 'Host'):
        host_taxonomy = []
        ncbi = NCBITaxa()
        for genus in rbp_embeddings[host_column].unique():
            taxonomy = [None, None, None, None, None, genus]
            try:
                genus_id = ncbi.get_name_translator([genus])[genus][0]
            except KeyError:
                try:
                    print(genus)
                    genus = 'candidatus ' + genus
                    genus_id = ncbi.get_name_translator([genus])[genus][0]
                except KeyError:
                    print(genus)
                    continue

            lineage_id = ncbi.get_lineage(genus_id)
            lineage = ncbi.get_rank(lineage_id)

            for taxon_id, rank in lineage.items():
                if rank == 'superkingdom':
                    taxonomy[0] = ncbi.get_taxid_translator([taxon_id])[taxon_id].lower()
                elif rank == 'phylum':
                    taxonomy[1] = ncbi.get_taxid_translator([taxon_id])[taxon_id].lower()
                elif rank == 'class':
                    taxonomy[2] = ncbi.get_taxid_translator([taxon_id])[taxon_id].lower()
                elif rank == 'order':
                    taxonomy[3] = ncbi.get_taxid_translator([taxon_id])[taxon_id].lower()
                elif rank == 'family':
                    taxonomy[4] = ncbi.get_taxid_translator([taxon_id])[taxon_id].lower()

            host_taxonomy.append(taxonomy)
            
        return host_taxonomy
    
    def get_sequences(self, rbps_with_accession, protein, *fasta_folders):
        util = RBPPredictionUtil()
        
        sequences = []
        for rbp in rbps_with_accession:
            entry = [rbp[0]]
            sequence = util.get_sequence(rbp[0], rbp[1], protein, *fasta_folders)
            
            if sequence is None:
                print(sequence)
                
            entry.append(sequence)
            sequences.append(entry)  
            
        return sequences
    
    # ===============
    # Embedding Files
    # ===============
    def get_rbp_embeddings(self, folder):
        rbps = []
            
        for file in os.listdir(f'{folder}'):
            with open(f'{folder}/{file}', 'r') as embeddings:
                for embedding in embeddings:
                    embedding_list = embedding.rstrip('\n').split(',')
                    rbp_id = embedding_list[0]
                    embedding_vals = embedding_list[1:]
                    if rbp_id != 'ID':
                        rbps.append([rbp_id] + embedding_vals)
                        
        return rbps
    
    def get_rbp_embeddings_df(self, plm, folder):
        if plm == 'PROTTRANSBERT':
            folder += '/master'
            
        rbp_df = pd.DataFrame(self.get_rbp_embeddings(folder))
        rbp_df.rename(columns = {0: 'Protein ID'}, inplace = True)
        
        rbp_df2 = rbp_df.loc[:, rbp_df.columns != 'Protein ID']
        rbp_df2 = rbp_df2.astype(np.float64)
        
        rbp_df2['Protein ID'] = rbp_df['Protein ID'].values
        col1 = rbp_df2.pop('Protein ID')
        rbp_df2.insert(0, 'Protein ID', col1)
        
        return rbp_df2
        
    # ==============
    # Classification
    # ==============
    def random_train_test_split(self, rbp_embeddings, taxon, embeddings_size = None, feature_columns = None):
        if not feature_columns:
            feature_columns = [str(i) for i in range(1, embeddings_size + 1)]
        
        X = rbp_embeddings.loc[:, rbp_embeddings.columns.isin(feature_columns)]
        y = rbp_embeddings.loc[:, rbp_embeddings.columns.isin([taxon])]

        X_train, X_test, y_train, y_test = train_test_split(X, y, stratify = y, test_size = 0.3, random_state = self.RANDOM_NUM)
        
        all_counts = rbp_embeddings[taxon].value_counts()
        pre_counts = y_train[taxon].value_counts()
        post_counts = y_test[taxon].value_counts()
        
        counts = []
        for rank in all_counts.index:
            entry = [rank]
            try:
                entry.append(pre_counts[rank])
            except KeyError:
                entry.append(0)

            try:
                entry.append(post_counts[rank])
            except KeyError:
                entry.append(0)

            try:
                entry.append(all_counts[rank])
            except KeyError:
                entry.append(0)

            counts.append(entry)
            
        return counts, X_train, X_test, y_train, y_test
    
    def get_unknown_hosts(self, rbp_embeddings, taxon, embeddings_size = None, feature_columns = None):
        if not feature_columns:
            feature_columns = [str(i) for i in range(1, embeddings_size + 1)]
        
        X = rbp_embeddings.loc[:, rbp_embeddings.columns.isin(feature_columns)]
        y = rbp_embeddings.loc[:, rbp_embeddings.columns.isin([taxon])]
        
        y[taxon] = y[taxon].apply(lambda x: 'unknown')

        return X, y
    
    def xy_split(self, rbp_embeddings, taxon, embeddings_size = None, feature_columns = None):
        if not feature_columns:
            feature_columns = [str(i) for i in range(1, embeddings_size + 1)]
        
        X = rbp_embeddings.loc[:, rbp_embeddings.columns.isin(feature_columns)]
        y = rbp_embeddings.loc[:, rbp_embeddings.columns.isin([taxon])]
        
        return X, y
    
    # ========
    # PhagesDB
    # ========
    
    def get_names_no_accession(self, start, end, phages, not_in_genbank_dir):
        phage_names = set()

        for i in range(start, end + 1):
            if str(i) in phages:
                for record in SeqIO.parse(f'{not_in_genbank_dir}/{str(i)}.ffn', 'fasta'):
                    phage_names.add(record.description.split(' ')[1])
                    
        return phage_names
    
    def get_names_no_accession_with_accession(self, start, end, phages, not_in_genbank_dir):
        phage_names = {}

        for i in range(start, end + 1):
            if str(i) in phages:
                for record in SeqIO.parse(f'{not_in_genbank_dir}/{str(i)}.ffn', 'fasta'):
                    phage_names[record.description.split(' ')[1]] = record.description.split(' ')[0]
                    
        return phage_names
    
    def predict_with_threshold(self, proba, y_test, y_pred, unknown_threshold = 0, display = False):
        y_pred_copy = copy.deepcopy(y_pred)

        unknown = 'unknown'
        
        largest = []
        second_largest = []
        
        for row in proba:
            row_sorted = sorted(row)
            largest.append(row_sorted[-1])
            second_largest.append(row_sorted[-2])
        
        for idx, (largest_val, second_largest_val) in enumerate(zip(largest, second_largest)):
            if largest_val - second_largest_val  < unknown_threshold:
                y_pred_copy[idx] = unknown

        if display:
            print('Unknown threshold:', str(unknown_threshold * 100) + '%')
            print(classification_report(y_test, y_pred_copy, digits = 4))
            print('===================')
            
        return (precision_recall_fscore_support(y_test, y_pred_copy, average = None),
                precision_recall_fscore_support(y_test, y_pred_copy, average = 'micro'),
                precision_recall_fscore_support(y_test, y_pred_copy, average = 'macro'),
                precision_recall_fscore_support(y_test, y_pred_copy, average = 'weighted'),
                proba,
                y_test,
                y_pred,
                y_pred_copy)