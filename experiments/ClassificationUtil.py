"""
==============================================================
This script contains the utility functions for classification.
@author    Mark Edward M. Gonzales
==============================================================
"""

import os
import copy
import math
import pickle

import pandas as pd
import numpy as np

from ete3 import NCBITaxa
from Bio import SeqIO, Align
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Align import substitution_matrices
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, precision_recall_fscore_support
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils.fixes import loguniform
from scipy.stats import randint
from joblib import dump
from collections import defaultdict

from RBPPredictionUtil import RBPPredictionUtil
from ConstantsUtil import ConstantsUtil
from SequenceParsingUtil import SequenceParsingUtil

class ClassificationUtil(object):
    def __init__(self, complete_embeddings_dir = None, RANDOM_NUM = 42):
        """
        Constructor
        
        Parameters:
        - complete_embeddings_dir: File path of the directory containing the embeddings of the RBPs
        - RANDOM_NUM: Random number for reproducibility
        """
        self.complete_embeddings_dir = complete_embeddings_dir
        self.RANDOM_NUM = RANDOM_NUM
        self.inphared_gb = None
    
    # =============
    # Preprocessing
    # =============
    def set_inphared_gb(self, inphared_gb):
        """
        Sets the consolidated GenBank entries of the entries fetched via INPHARED
        
        Parameters:
        - inphared_gb: File path of the consolidated GenBank entries of the entries fetched via INPHARED
        """
        self.inphared_gb = inphared_gb
    
    def get_phages(self):
        """
        Retrieves the phage IDs of the phages in the dataset
        
        Returns:
        - Set of phage IDs of the phages in the dataset
        """
        phages = set()
        for file in os.listdir(f'{self.complete_embeddings_dir}'):
            phages.add(file[:-len('-rbp-embeddings.csv')])
            
        return phages
    
    def get_rbps(self):
        """
        Retrieves the protein IDs of the RBPs in the dataset
        
        Returns:
        - Set of protein IDs of the RBPs in the dataset
        """
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
        """
        Retrieves the taxonomical information of the hosts (particularly the superkingdom, phylum, class, order, and family)
        since INPHARED returns only the hosts at genus level
        
        Parameters:
        - rbp_embeddings: DataFrame containing the RBP embeddings
        - host_column: Name of the host genus
        
        Returns:
        - List containing the taxonomical information of the hosts; each item in the list corresponds to one host
        """
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
        """
        Retrieves the RBP protein (or nucleotide) sequences
        
        Parameters:
        - rbps_with_accession: DataFrame containing the RBP embeddings
        - protein: Protein ID
        - fasta_folders: File paths of the directories containing the FASTA (or FFN) files with the sequences
        
        Returns:
        - List containing the RBP protein (or nucleotide) sequences; each item in the list corresponds to one RBP 
        """
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
        """
        Retrieves the protein IDs and embeddings of the RBPs in the given directory
        
        Parameters:
        - folder: File path of the directory containing the embeddings of the RBPs
        
        Returns:
        - List containing the protein IDs and embeddings of the RBPs in the directory; each item in the list corresponds to one RBP
        """
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
        """
        Constructs a DataFrame containing the protein IDs and embeddings of the RBPs in the given directory
        
        Parameters:
        - plm: Protein language model used to generate the embeddings
        - folder: File path of the directory containing the embeddings of the RBPs
        
        Returns:
        - DataFrame containing the protein IDs and embeddings of the RBPs in the directory
        """
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
        """
        Constructs the training and test sets for model training and evaluation
        
        Parameters:
        - rbp_embeddings: DataFrame containing the phages, hosts, and vector representations of the phages
        - taxon: Taxonomical level of classification
        - embeddings_size: Dimension of the embedding space
        - feature_columns: List of the column headers corresponding to the features
        
        Returns:
        - List containing the number of training and test samples for each class label; each item in the list corresponds to one class label
        - List containing the feature representations of the samples in the training set
        - List containing the class labels of the samples in the training set
        - List containing the feature representations of the samples in the test set
        - List containing the class labels of the samples in the test set
        """
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
        """
        Isolates the samples whose host is outside the class labels (i.e., outside the top 25% hosts)
        
        Parameters:
        - rbp_embeddings: DataFrame containing the phages, hosts, and vector representations of the phages
        - taxon: Taxonomical level of classification
        - embeddings_size: Dimension of the embedding space
        - feature_columns: List of the column headers corresponding to the features
        
        Returns:
        - List containing the feature representations of the samples whose host is outside the class labels 
          (i.e., outside the top 25% hosts)
        - List containing the class labels of the samples whose host is outside the class labels
          (i.e., outside the top 25% hosts). The class label of each of these samples is set to 'others'.
        """
        if not feature_columns:
            feature_columns = [str(i) for i in range(1, embeddings_size + 1)]
        
        X = rbp_embeddings.loc[:, rbp_embeddings.columns.isin(feature_columns)]
        y = rbp_embeddings.loc[:, rbp_embeddings.columns.isin([taxon])]
        
        y[taxon] = y[taxon].apply(lambda x: 'others')

        return X, y
    
    def xy_split(self, rbp_embeddings, taxon, embeddings_size = None, feature_columns = None):
        """
        Separates the feature representations and the class labels of the samples
        
        Parameters:
        - rbp_embeddings: DataFrame containing the phages, hosts, and vector representations of the phages
        - taxon: Taxonomical level of classification
        - embeddings_size: Dimension of the embedding space
        - feature_columns: List of the column headers corresponding to the features
        
        Returns:
        - List containing the feature representations of the samples
        - List containing the class labels of the samples
        """
        if not feature_columns:
            feature_columns = [str(i) for i in range(1, embeddings_size + 1)]
        
        X = rbp_embeddings.loc[:, rbp_embeddings.columns.isin(feature_columns)]
        y = rbp_embeddings.loc[:, rbp_embeddings.columns.isin([taxon])]
        
        return X, y
    
    def predict_with_threshold(self, proba, y_test, y_pred, unknown_threshold = 0, display = False):
        """
        Adjusts the predicted class labels in view of the given confidence threshold.
        In particular, a sample is classified under its predicted class label if and only if the difference between the largest
           and second largest class probabilities is greater than or equal to the given confidence threshold. Otherwise,
           the sample is classified as 'others' (i.e., falling outside the class labels).
        
        Parameters:
        - proba: Class probabilities
        - y_test: True class labels
        - y_pred: Predicted class labels
        - unknown_threshold: Confidence threshold
        - display: True if the per-class evaluation results are to be displayed; False, otherwise
        
        Returns:
        - List containing the per-class precision, recall, and F1
        - List containing the micro-precision, recall, and F1
        - List containing the macro-precision, recall, and F1
        - List containing the weighted precision, recall, and F1
        - Class probabilities
        - True class labels
        - Original predicted class labels 
        - Predicted class labels in view of the confidence threshold
        """
        y_pred_copy = copy.deepcopy(y_pred)

        unknown = 'others'
        
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
            print('Confidence threshold k:', str(unknown_threshold * 100) + '%')
            print(classification_report(y_test, y_pred_copy, digits = 4, zero_division = 0))
            print('===================')
            
        return (precision_recall_fscore_support(y_test, y_pred_copy, average = None, zero_division = 0),
                precision_recall_fscore_support(y_test, y_pred_copy, average = 'micro', zero_division = 0),
                precision_recall_fscore_support(y_test, y_pred_copy, average = 'macro', zero_division = 0),
                precision_recall_fscore_support(y_test, y_pred_copy, average = 'weighted', zero_division = 0),
                proba,
                y_test,
                y_pred,
                y_pred_copy)
    
    def classify(self, plm, save_feature_importance = False, display_feature_importance = False, feature_columns = None):
        """
        Trains a random forest model for phage-host interaction prediction and evaluates the model performance
        
        Parameters:
        - plm: Protein language model used to generate the embeddings
        - save_feature_importance: True if the components with the highest Gini importance are to be saved in a Pickle file;
                                   False, otherwise
        - display_feature_importance: True if the components with the highest Gini importance are to be displayed;
                                      False, otherwise
        - feature_columns: List of the column headers corresponding to the features
        
        Returns:
        - If display_feature_importance is set to True, the features with the highest Gini importance are returned
        - Otherwise, the function returns None
        """
        constants = ConstantsUtil()

        # Load data
        rbp_embeddings = pd.read_csv(f'{constants.INPHARED}/{constants.DATA}/{constants.PLM_EMBEDDINGS_CSV[plm]}', 
                                     low_memory = False)
        rbp_embeddings['Modification Date'] = pd.to_datetime(rbp_embeddings['Modification Date'])

        # Get only the top 25% hosts
        all_counts = rbp_embeddings['Host'].value_counts()
        TOP_X_PERCENT = 0.25
        top_x = math.floor(all_counts.shape[0] * TOP_X_PERCENT)

        top_genus = set()
        genus_counts = all_counts.index
        for entry in genus_counts[:top_x]:
            top_genus.add(entry)

        # Construct the training and test sets
        print("Constructing training and test sets...")

        rbp_embeddings_top = rbp_embeddings[rbp_embeddings['Host'].isin(top_genus)]

        counts = X_train = X_test = y_train = y_test = None
        if feature_columns is None:
            counts, X_train, X_test, y_train, y_test = self.random_train_test_split(rbp_embeddings_top, 'Host',
                                                                                    embeddings_size = rbp_embeddings.shape[1] - constants.INPHARED_EXTRA_COLS)
        else:
            counts, X_train, X_test, y_train, y_test = self.random_train_test_split(rbp_embeddings_top, 'Host',
                                                                                    feature_columns = feature_columns)

        counts_df = pd.DataFrame(counts, columns = ['Genus', f'Train', f'Test', 'Total'])

        unknown_hosts_X = unknown_hosts_y = None
        if feature_columns is None:
            unknown_hosts_X, unknown_hosts_y = self.get_unknown_hosts(rbp_embeddings[~rbp_embeddings['Host'].isin(top_genus)], 'Host',
                                                                      embeddings_size = rbp_embeddings.shape[1] - constants.INPHARED_EXTRA_COLS)
        else:
            unknown_hosts_X, unknown_hosts_y = self.get_unknown_hosts(rbp_embeddings[~rbp_embeddings['Host'].isin(top_genus)], 'Host',
                                                                      feature_columns = feature_columns)

        X_test = X_test.append(unknown_hosts_X)
        y_test = y_test.append(unknown_hosts_y)

        # Construct the model
        print("Training the model...")

        clf = RandomForestClassifier(random_state = self.RANDOM_NUM, class_weight = 'balanced',
                                     max_features = 'sqrt',
                                     min_samples_leaf = 1,
                                     min_samples_split = 2,
                                     n_estimators = 150,
                                     n_jobs = -1)

        clf.fit(X_train, y_train.values.ravel())
        y_pred = clf.predict(X_test)
        proba = clf.predict_proba(X_test)

        # Save the results
        print("Saving evaluation results...")

        results = []
        for threshold in range(0, 101, 10):
            results.append(self.predict_with_threshold(proba, y_test, y_pred, unknown_threshold = threshold / 100, display = True))

        if not os.path.exists(constants.TEMP_RESULTS):
            os.makedirs(constants.TEMP_RESULTS)

        with open(constants.PLM_RESULTS[plm], 'wb') as f:
            pickle.dump(results, f)

        # Save the trained model
        if not os.path.exists(constants.TRAINED_MODEL):
            os.makedirs(constants.TRAINED_MODEL)

        dump(clf, constants.PLM_TRAINED_MODEL[plm])

        # Only for the best-performing model (ProtT5)
        if save_feature_importance:
            # Save feature importance
            feature_importances = clf.feature_importances_
            # Add 1 so that component 1 is mapped to 1 (instead of 0)
            indices = [x + 1 for _, x in sorted(zip(feature_importances, range(len(feature_importances))), reverse = True)]
                
            with open(constants.FEATURE_IMPORTANCE, 'wb') as f:
                pickle.dump(indices, f)
                
        if display_feature_importance:
            feature_importances = clf.feature_importances_
            indices = [x for _, x in sorted(zip(feature_importances, range(len(feature_importances))), reverse = True)]
            names = [feature_columns[index] for index in indices]
            print(names)
            
            return names
            
        # Display progress
        print("Finished")
        print("==============")
        
    def classify_handpicked_embeddings(self, feature_columns, filename):
        """
        Trains a random forest model for phage-host interaction prediction and evaluates the model performance.
        This method is for investigating the change in performance when the ProtT5 embeddings are combined with handcrafted features.
        
        Parameters:
        - feature_columns: List of the column headers corresponding to the features
        - filename: File name of the Pickle file in which the evaluation results are to be saved
        """
        constants = ConstantsUtil()
        
        # Load data
        boeckaerts = pd.read_csv(f'{constants.INPHARED}/{constants.DATA}/{constants.PLM_EMBEDDINGS_CSV["BOECKAERTS"]}',
                                 low_memory = False)
        prott5 = pd.read_csv(f'{constants.INPHARED}/{constants.DATA}/{constants.PLM_EMBEDDINGS_CSV["PROTT5"]}',
                             low_memory = False)
        rbp_embeddings = pd.merge(boeckaerts, prott5, how = 'inner', validate = 'one_to_one')
        
        # Get only the top 25% hosts
        all_counts = rbp_embeddings['Host'].value_counts()
        TOP_X_PERCENT = 0.25
        top_x = math.floor(all_counts.shape[0] * TOP_X_PERCENT)

        top_genus = set()
        genus_counts = all_counts.index
        for entry in genus_counts[:top_x]:
            top_genus.add(entry)
            
        merged_top = rbp_embeddings[rbp_embeddings['Host'].isin(top_genus)]
        
        # Construct the training and test sets
        print("Constructing training and test sets...")

        counts, X_train, X_test, y_train, y_test = self.random_train_test_split(merged_top, 'Host',
                                                                                feature_columns = feature_columns + [str(i) for i in range(1, prott5.shape[1] - constants.INPHARED_EXTRA_COLS + 1)])

        unknown_hosts_X, unknown_hosts_y = self.get_unknown_hosts(rbp_embeddings[~rbp_embeddings['Host'].isin(top_genus)], 'Host',
                                                                  feature_columns = feature_columns + [str(i) for i in range(1, prott5.shape[1] - constants.INPHARED_EXTRA_COLS + 1)])
        X_test = X_test.append(unknown_hosts_X)
        y_test = y_test.append(unknown_hosts_y)
        
        # Construct the model
        print("Training the model...")

        clf = RandomForestClassifier(random_state = self.RANDOM_NUM, class_weight = 'balanced',
                                     max_features = 'sqrt',
                                     min_samples_leaf = 1,
                                     min_samples_split = 2,
                                     n_estimators = 150,
                                     n_jobs = -1)

        clf.fit(X_train, y_train.values.ravel())
        y_pred = clf.predict(X_test)
        proba = clf.predict_proba(X_test)

        # Save the results
        print("Saving evaluation results...")

        results = []
        for threshold in range(0, 101, 10):
            results.append(self.predict_with_threshold(proba, y_test, y_pred, unknown_threshold = threshold / 100, display = False))

        if not os.path.exists(constants.TEMP_RESULTS):
            os.makedirs(constants.TEMP_RESULTS)

        with open(f'{constants.TEMP_RESULTS}/prott5_{filename}.pickle', 'wb') as f:
            pickle.dump(results, f)
    
    # =============
    # Miscellaneous
    # =============
    def get_train_test_sets(self):
        """
        Performs a stratified train-test split. The train-test partition returned by this method is the same as the one used
        in the methods for classification.
        
        Returns:
        - Training set (features)
        - Training set (labels / host genera)
        - Test set (features)
        - Test set (labels / host genera)
        """
        constants = ConstantsUtil()
        plm = 'PROTT5'

        # Load data
        rbp_embeddings = pd.read_csv(f'{constants.INPHARED}/{constants.DATA}/{constants.PLM_EMBEDDINGS_CSV[plm]}', 
                                     low_memory = False)
        rbp_embeddings['Modification Date'] = pd.to_datetime(rbp_embeddings['Modification Date'])

        # Get only the top 25% hosts
        all_counts = rbp_embeddings['Host'].value_counts()
        TOP_X_PERCENT = 0.25
        top_x = math.floor(all_counts.shape[0] * TOP_X_PERCENT)

        top_genus = set()
        genus_counts = all_counts.index
        for entry in genus_counts[:top_x]:
            top_genus.add(entry)

        # Construct the training and test sets
        print("Constructing training and test sets...")

        rbp_embeddings_top = rbp_embeddings[rbp_embeddings['Host'].isin(top_genus)]

        counts, X_train, X_test, y_train, y_test = self.random_train_test_split(rbp_embeddings_top, 'Host',
                                                                                embeddings_size = rbp_embeddings.shape[1] - constants.INPHARED_EXTRA_COLS)

        counts_df = pd.DataFrame(counts, columns = ['Genus', f'Train', f'Test', 'Total'])

        unknown_hosts_X, unknown_hosts_y = self.get_unknown_hosts(rbp_embeddings[~rbp_embeddings['Host'].isin(top_genus)], 'Host',
                                                                  embeddings_size = rbp_embeddings.shape[1] - constants.INPHARED_EXTRA_COLS)

        X_test = X_test.append(unknown_hosts_X)
        y_test = y_test.append(unknown_hosts_y)
        
        return X_train, y_train, X_test, y_test