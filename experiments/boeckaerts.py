"""
COLLECTION OF FUNCTIONS FOR PROTEIN SEQUENCE FEATURE CONSTRUCTION & BLAST PREDICTION

Created on Thu Nov  9 13:29:44 2017

@author: dimiboeckaerts

Some of the code below is taken from the following Github repo:
https://github.com/Superzchen/iFeature
(Chen et al., 2018. Bioinformatics.)
"""

# IMPORT LIBRARIES
# --------------------------------------------------
import os
import math
import warnings
import numpy as np
import scipy as sp
import datetime as dt
from numba import jit
import matplotlib.pyplot as plt
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Entrez, pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from sklearn.preprocessing import label_binarize
from sklearn.exceptions import UndefinedMetricWarning
from Bio.Blast.Applications import NcbiblastpCommandline
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.metrics import accuracy_score, make_scorer, f1_score, precision_score, recall_score
from sklearn.metrics import precision_recall_curve, auc


# DNA FEATURES
# --------------------------------------------------
def dna_features(dna_sequences):
    """
    This function calculates a variety of properties from a DNA sequence.
    
    Input: a list of DNA sequence (can also be length of 1)
    Output: a dataframe of features
    """
    
    import numpy as np
    import pandas as pd
    from Bio.SeqUtils import GC, CodonUsage

    A_freq = []; T_freq = []; C_freq = []; G_freq = []; GC_content = []
    codontable = {'ATA':[], 'ATC':[], 'ATT':[], 'ATG':[], 'ACA':[], 'ACC':[], 'ACG':[], 'ACT':[],
    'AAC':[], 'AAT':[], 'AAA':[], 'AAG':[], 'AGC':[], 'AGT':[], 'AGA':[], 'AGG':[],
    'CTA':[], 'CTC':[], 'CTG':[], 'CTT':[], 'CCA':[], 'CCC':[], 'CCG':[], 'CCT':[],
    'CAC':[], 'CAT':[], 'CAA':[], 'CAG':[], 'CGA':[], 'CGC':[], 'CGG':[], 'CGT':[],
    'GTA':[], 'GTC':[], 'GTG':[], 'GTT':[], 'GCA':[], 'GCC':[], 'GCG':[], 'GCT':[],
    'GAC':[], 'GAT':[], 'GAA':[], 'GAG':[], 'GGA':[], 'GGC':[], 'GGG':[], 'GGT':[],
    'TCA':[], 'TCC':[], 'TCG':[], 'TCT':[], 'TTC':[], 'TTT':[], 'TTA':[], 'TTG':[],
    'TAC':[], 'TAT':[], 'TAA':[], 'TAG':[], 'TGC':[], 'TGT':[], 'TGA':[], 'TGG':[]}
    
    for item in dna_sequences:
        # nucleotide frequencies
        A_freq.append(item.count('A')/len(item))
        T_freq.append(item.count('T')/len(item))
        C_freq.append(item.count('C')/len(item))
        G_freq.append(item.count('G')/len(item))
    
        # GC content
        GC_content.append(GC(item))
    
        # codon frequency: count codons, normalize counts, add to dict
        codons = [item[i:i+3] for i in range(0, len(item), 3)]
        l = []
        for key in codontable.keys():
            l.append(codons.count(key))
        l_norm = [float(i)/sum(l) for i in l]
        
        for j, key in enumerate(codontable.keys()):
            codontable[key].append(l_norm[j])
     
    # codon usage bias (_b)
    synonym_codons = CodonUsage.SynonymousCodons
    codontable2 = {'ATA_b':[], 'ATC_b':[], 'ATT_b':[], 'ATG_b':[], 'ACA_b':[], 'ACC_b':[], 'ACG_b':[], 'ACT_b':[],
    'AAC_b':[], 'AAT_b':[], 'AAA_b':[], 'AAG_b':[], 'AGC_b':[], 'AGT_b':[], 'AGA_b':[], 'AGG_b':[],
    'CTA_b':[], 'CTC_b':[], 'CTG_b':[], 'CTT_b':[], 'CCA_b':[], 'CCC_b':[], 'CCG_b':[], 'CCT_b':[],
    'CAC_b':[], 'CAT_b':[], 'CAA_b':[], 'CAG_b':[], 'CGA_b':[], 'CGC_b':[], 'CGG_b':[], 'CGT_b':[],
    'GTA_b':[], 'GTC_b':[], 'GTG_b':[], 'GTT_b':[], 'GCA_b':[], 'GCC_b':[], 'GCG_b':[], 'GCT_b':[],
    'GAC_b':[], 'GAT_b':[], 'GAA_b':[], 'GAG_b':[], 'GGA_b':[], 'GGC_b':[], 'GGG_b':[], 'GGT_b':[],
    'TCA_b':[], 'TCC_b':[], 'TCG_b':[], 'TCT_b':[], 'TTC_b':[], 'TTT_b':[], 'TTA_b':[], 'TTG_b':[],
    'TAC_b':[], 'TAT_b':[], 'TAA_b':[], 'TAG_b':[], 'TGC_b':[], 'TGT_b':[], 'TGA_b':[], 'TGG_b':[]}

    for item1 in dna_sequences:
        codons = [item1[l:l+3] for l in range(0, len(item1), 3)]
        codon_counts = []
    
        # count codons corresponding to codontable (not codontable2 because keynames changed!)
        for key in codontable.keys():
            codon_counts.append(codons.count(key))
        
        # count total for synonymous codons, divide each synonym codon count by total
        for key_syn in synonym_codons.keys():
            total = 0
            for item2 in synonym_codons[key_syn]:
                total += codons.count(item2)
            for j, key_table in enumerate(codontable.keys()):
                if (key_table in synonym_codons[key_syn]) & (total != 0):
                    codon_counts[j] /= total
                
        # add corrected counts to codontable2 (also corresponds to codontable which was used to count codons)
        for k, key_table in enumerate(codontable2.keys()):
            codontable2[key_table].append(codon_counts[k])
            
    # make new dataframes & standardize
    features_codonbias = pd.DataFrame.from_dict(codontable2)
    features_dna = pd.DataFrame.from_dict(codontable)
    features_dna['A_freq'] = np.asarray(A_freq)
    features_dna['T_freq'] = np.asarray(T_freq)
    features_dna['C_freq'] = np.asarray(C_freq)
    features_dna['G_freq'] = np.asarray(G_freq)
    features_dna['GC'] = np.asarray(GC_content)
    
    # concatenate dataframes & return
    features = pd.concat([features_dna, features_codonbias], axis=1)
    return features


# PROTEIN FEATURE: BASICS
# --------------------------------------------------
def protein_features(protein_sequences):
    """
    This function calculates a number of basic properties for a list of protein sequences
    
    Input: list of protein sequences (as strings), length can also be 1
    Output: a dataframe of features
    """
    
    import numpy as np
    import pandas as pd
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    
    # AA frequency and protein characteristics
    mol_weight = []; aromaticity = []; instability = []; flexibility = []; prot_length = []
    pI = []; helix_frac = []; turn_frac = []; sheet_frac = []
    frac_aliph = []; frac_unch_polar = []; frac_polar = []; frac_hydrophob = []; frac_pos = []; frac_sulfur = []
    frac_neg = []; frac_amide = []; frac_alcohol = []
    AA_dict = {'G': [], 'A': [], 'V': [], 'L': [], 'I': [], 'F': [], 'P': [], 'S': [], 'T': [], 'Y': [],
           'Q': [], 'N': [], 'E': [], 'D': [], 'W': [], 'H': [], 'R': [], 'K': [], 'M': [], 'C': []}
    
    for item in protein_sequences:
        # calculate various protein properties
        prot_length.append(len(item))
        frac_aliph.append((item.count('A')+item.count('G')+item.count('I')+item.count('L')+item.count('P')
                       +item.count('V'))/len(item))
        frac_unch_polar.append((item.count('S')+item.count('T')+item.count('N')+item.count('Q'))/len(item))
        frac_polar.append((item.count('Q')+item.count('N')+item.count('H')+item.count('S')+item.count('T')+item.count('Y')
                      +item.count('C')+item.count('M')+item.count('W'))/len(item))
        frac_hydrophob.append((item.count('A')+item.count('G')+item.count('I')+item.count('L')+item.count('P')
                        +item.count('V')+item.count('F'))/len(item))
        frac_pos.append((item.count('H')+item.count('K')+item.count('R'))/len(item))
        frac_sulfur.append((item.count('C')+item.count('M'))/len(item))
        frac_neg.append((item.count('D')+item.count('E'))/len(item))
        frac_amide.append((item.count('N')+item.count('Q'))/len(item))
        frac_alcohol.append((item.count('S')+item.count('T'))/len(item))
        protein_chars = ProteinAnalysis(item) 
        mol_weight.append(protein_chars.molecular_weight())
        aromaticity.append(protein_chars.aromaticity())
        instability.append(protein_chars.instability_index())
        flexibility.append(np.mean(protein_chars.flexibility()))
        pI.append(protein_chars.isoelectric_point())
        H, T, S = protein_chars.secondary_structure_fraction()
        helix_frac.append(H)
        turn_frac.append(T)
        sheet_frac.append(S)
    
        # calculate AA frequency
        for key in AA_dict.keys():
            AA_dict[key].append(item.count(key)/len(item))
            
    # make new dataframe & return
    features_protein = pd.DataFrame.from_dict(AA_dict)
    features_protein['protein_length'] = np.asarray(prot_length)
    features_protein['mol_weight'] = np.asarray(mol_weight)
    features_protein['aromaticity'] = np.asarray(aromaticity)
    features_protein['instability'] = np.asarray(instability)
    features_protein['flexibility'] = np.asarray(flexibility)
    features_protein['pI'] = np.asarray(pI)
    features_protein['frac_aliphatic'] = np.asarray(frac_aliph)
    features_protein['frac_uncharged_polar'] = np.asarray(frac_unch_polar)
    features_protein['frac_polar'] = np.asarray(frac_polar)
    features_protein['frac_hydrophobic'] = np.asarray(frac_hydrophob)
    features_protein['frac_positive'] = np.asarray(frac_pos)
    features_protein['frac_sulfur'] = np.asarray(frac_sulfur)
    features_protein['frac_negative'] = np.asarray(frac_neg)
    features_protein['frac_amide'] = np.asarray(frac_amide)
    features_protein['frac_alcohol'] = np.asarray(frac_alcohol)
    features_protein['AA_frac_helix'] = np.asarray(helix_frac)
    features_protein['AA_frac_turn'] = np.asarray(turn_frac)
    features_protein['AA_frac_sheet'] = np.asarray(sheet_frac)
    
    return features_protein


# PROTEIN FEATURE: COMPOSITION
# --------------------------------------------------
def Count1(seq1, seq2):
	sum = 0
	for aa in seq1:
		sum = sum + seq2.count(aa)
	return sum

def CTDC(sequence):
    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'hydrophobicity_ARGP820101': 'QSTNGDE',
        'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
        'hydrophobicity_PONP930101': 'KPDESNQT',
        'hydrophobicity_CASG920101': 'KDEQPSRNTG',
        'hydrophobicity_ENGD860101': 'RDKENQHYP',
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity':        'LIFWCMVY',
        'polarizability':  'GASDT',
        'charge':          'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess':   'ALFCGIVW'
    }
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'hydrophobicity_ARGP820101': 'RAHCKMV',
        'hydrophobicity_ZIMJ680101': 'HMCKV',
        'hydrophobicity_PONP930101': 'GRHA',
        'hydrophobicity_CASG920101': 'AHYMLV',
        'hydrophobicity_ENGD860101': 'SGTAW',
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity':        'PATGS',
        'polarizability':  'CPNVEQIL',
        'charge':          'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess':   'RKQEND'
    }

    property = ['solventaccess'] # can be expanded to include other properties
    encoding = []

    for p in property:
        c1 = Count1(group1[p], sequence) / len(sequence)
        c2 = Count1(group2[p], sequence) / len(sequence)
        c3 = 1 - c1 - c2
        encoding = encoding + [c1, c2, c3]
        
    return encoding


# PROTEIN FEATURE: TRANSITION
# --------------------------------------------------
def CTDT(sequence):
    """
    Every number in the encoding tells us as a percentage over all AA pairs, how many transitions
    occured from group x to y or vice versa for a specific property p.
    """
    group1 = {'hydrophobicity_PRAM900101': 'RKEDQN','hydrophobicity_ARGP820101': 'QSTNGDE',
              'hydrophobicity_ZIMJ680101': 'QNGSWTDERA','hydrophobicity_PONP930101': 'KPDESNQT',
              'hydrophobicity_CASG920101': 'KDEQPSRNTG','hydrophobicity_ENGD860101': 'RDKENQHYP',
              'hydrophobicity_FASG890101': 'KERSQD','normwaalsvolume': 'GASTPDC','polarity': 'LIFWCMVY',
              'polarizability': 'GASDT','charge':'KR', 'secondarystruct': 'EALMQKRH','solventaccess': 'ALFCGIVW'}
    
    group2 = {'hydrophobicity_PRAM900101': 'GASTPHY','hydrophobicity_ARGP820101': 'RAHCKMV',
           'hydrophobicity_ZIMJ680101': 'HMCKV','hydrophobicity_PONP930101': 'GRHA',
           'hydrophobicity_CASG920101': 'AHYMLV','hydrophobicity_ENGD860101': 'SGTAW',
           'hydrophobicity_FASG890101': 'NTPG','normwaalsvolume': 'NVEQIL','polarity': 'PATGS',
           'polarizability': 'CPNVEQIL','charge': 'ANCQGHILMFPSTWYV', 
           'secondarystruct': 'VIYCWFT', 'solventaccess': 'RKQEND'}
    
    group3 = {'hydrophobicity_PRAM900101': 'CLVIMFW','hydrophobicity_ARGP820101': 'LYPFIW',
           'hydrophobicity_ZIMJ680101': 'LPFYI','hydrophobicity_PONP930101': 'YMFWLCVI',
           'hydrophobicity_CASG920101': 'FIWC','hydrophobicity_ENGD860101': 'CVLIMF',
           'hydrophobicity_FASG890101': 'AYHWVMFLIC','normwaalsvolume': 'MHKFRYW',
           'polarity': 'HQRKNED','polarizability': 'KMHFRYW','charge': 'DE',
           'secondarystruct': 'GNPSD','solventaccess': 'MSPTHY'}
    
    property = ['hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess']
    
    encoding = []
    aaPair = [sequence[j:j+2] for j in range(len(sequence)-1)]
    
    for p in property:
        c1221, c1331, c2332 = 0, 0, 0
        for pair in aaPair:
            if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                c1221 += 1
                continue
            if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
                c1331 += 1
                continue
            if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                c2332 += 1
        encoding.append(c1221/len(aaPair))
        encoding.append(c1331/len(aaPair))
        encoding.append(c2332/len(aaPair))
    
    return encoding


# PROTEIN FEATURE: DISTRIBUTION
# --------------------------------------------------
def Count2(aaSet, sequence):
	number = 0
	for aa in sequence:
		if aa in aaSet:
			number = number + 1
	cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
	cutoffNums = [i if i >=1 else 1 for i in cutoffNums]

	code = []
	for cutoff in cutoffNums:
		myCount = 0
		for i in range(len(sequence)):
			if sequence[i] in aaSet:
				myCount += 1
				if myCount == cutoff:
					code.append((i + 1) / len(sequence) * 100)
					break
		if myCount == 0:
			code.append(0)
	return code

def CTDD(sequence):
	group1 = {
		'hydrophobicity_PRAM900101': 'RKEDQN',
		'hydrophobicity_ARGP820101': 'QSTNGDE',
		'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
		'hydrophobicity_PONP930101': 'KPDESNQT',
		'hydrophobicity_CASG920101': 'KDEQPSRNTG',
		'hydrophobicity_ENGD860101': 'RDKENQHYP',
		'hydrophobicity_FASG890101': 'KERSQD',
		'normwaalsvolume': 'GASTPDC',
		'polarity':        'LIFWCMVY',
		'polarizability':  'GASDT',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	}
	group2 = {
		'hydrophobicity_PRAM900101': 'GASTPHY',
		'hydrophobicity_ARGP820101': 'RAHCKMV',
		'hydrophobicity_ZIMJ680101': 'HMCKV',
		'hydrophobicity_PONP930101': 'GRHA',
		'hydrophobicity_CASG920101': 'AHYMLV',
		'hydrophobicity_ENGD860101': 'SGTAW',
		'hydrophobicity_FASG890101': 'NTPG',
		'normwaalsvolume': 'NVEQIL',
		'polarity':        'PATGS',
		'polarizability':  'CPNVEQIL',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
	}
	group3 = {
		'hydrophobicity_PRAM900101': 'CLVIMFW',
		'hydrophobicity_ARGP820101': 'LYPFIW',
		'hydrophobicity_ZIMJ680101': 'LPFYI',
		'hydrophobicity_PONP930101': 'YMFWLCVI',
		'hydrophobicity_CASG920101': 'FIWC',
		'hydrophobicity_ENGD860101': 'CVLIMF',
		'hydrophobicity_FASG890101': 'AYHWVMFLIC',
		'normwaalsvolume': 'MHKFRYW',
		'polarity':        'HQRKNED',
		'polarizability':  'KMHFRYW',
		'charge':          'DE',
		'secondarystruct': 'GNPSD',
		'solventaccess':   'MSPTHY'
	}

	property = ['hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess']

	for p in property:
		encoding = Count2(group1[p], sequence) + Count2(group2[p], sequence) + Count2(group3[p], sequence)
	return encoding


# PROTEIN FEATURE: Z-SCALE
# --------------------------------------------------
def zscale(sequence):
    zdict = {
		'A': [0.24,  -2.32,  0.60, -0.14,  1.30], # A
		'C': [0.84,  -1.67,  3.71,  0.18, -2.65], # C
		'D': [3.98,   0.93,  1.93, -2.46,  0.75], # D
		'E': [3.11,   0.26, -0.11, -0.34, -0.25], # E
		'F': [-4.22,  1.94,  1.06,  0.54, -0.62], # F
		'G': [2.05,  -4.06,  0.36, -0.82, -0.38], # G
		'H': [2.47,   1.95,  0.26,  3.90,  0.09], # H
		'I': [-3.89, -1.73, -1.71, -0.84,  0.26], # I
		'K': [2.29,   0.89, -2.49,  1.49,  0.31], # K
		'L': [-4.28, -1.30, -1.49, -0.72,  0.84], # L
		'M': [-2.85, -0.22,  0.47,  1.94, -0.98], # M
		'N': [3.05,   1.62,  1.04, -1.15,  1.61], # N
		'P': [-1.66,  0.27,  1.84,  0.70,  2.00], # P
		'Q': [1.75,   0.50, -1.44, -1.34,  0.66], # Q
		'R': [3.52,   2.50, -3.50,  1.99, -0.17], # R
		'S': [2.39,  -1.07,  1.15, -1.39,  0.67], # S
		'T': [0.75,  -2.18, -1.12, -1.46, -0.40], # T
		'V': [-2.59, -2.64, -1.54, -0.85, -0.02], # V
		'W': [-4.36,  3.94,  0.59,  3.44, -1.59], # W
		'Y': [-2.54,  2.44,  0.43,  0.04, -1.47], # Y
		'-': [0.00,   0.00,  0.00,  0.00,  0.00], # -
	}
    
    z1, z2, z3, z4, z5 = 0, 0, 0, 0, 0
    for aa in sequence:
        z1 += zdict[aa][0]
        z2 += zdict[aa][1]
        z3 += zdict[aa][2]
        z4 += zdict[aa][3]
        z5 += zdict[aa][4]
    encoding = [z1/len(sequence), z2/len(sequence), z3/len(sequence), z4/len(sequence), z5/len(sequence)]
    
    return encoding


# GROUPING FOR GROUP K FOLD CV
# --------------------------------------------------
def define_groups(matrix, threshold):
    # make groups of sequences under identity threshold
    groups_dict = {}
    groups_list = []
    for i in range(matrix.shape[0]):
        for j in range(i+1, matrix.shape[1]):
            if (matrix[i,j] >= threshold) & (i not in groups_list) & (j not in groups_list) & (i not in groups_dict.keys()):
                groups_dict[i] = [j]
                groups_list.append(j)
            elif (matrix[i,j] >= threshold) & (i not in groups_list) & (j not in groups_list) & (i in groups_dict.keys()):
                groups_dict[i].append(j)
                groups_list.append(j)
                
    # MARK
    print(groups_dict)
    print(groups_list)

    # assign group numbers to clusters
    groups_array = np.zeros(matrix.shape[0])
    groups_index = 1
    for key in groups_dict.keys():
        groups_array[key] = groups_index
        for item in groups_dict[key]:
            groups_array[item] = groups_index
        groups_index += 1

    # assign group numbers to leftover sequences not in any cluster
    for i,item in enumerate(groups_array):
        if item == 0:
            groups_array[i] = groups_index
            groups_index += 1
    
    return groups_array


# NESTED CROSS-VALIDATION WITH GROUPKFOLD
# --------------------------------------------------
def NestedGroupKFold(model, X, y, parameter_grid, groups, class_weights, scorer=make_scorer(accuracy_score), 
                     inner_cv=GroupKFold(n_splits=4), outer_cv=GroupKFold(n_splits=4)):
    """
    Implements a nested version of GroupKFold cross-validation using GridSearchCV to evaluate models 
    that need hyperparameter tuning in settings where different groups exist in the available data.
    
    Dependencies: sklearn.model_selection, sklearn.metrics, numpy
    
    Input:
    - X, y: features and labels (must be NumPy arrays).
    - model, parameter_grid: the model instance and its parameter grid to be optimized.
    - groups: the groups to use in both inner- and outer loop.
    - scorer: the scoring to use in inner loop (default: accuracy).
    - inner_cv, outer_cv: the iterators for both CV-loops (default: GroupKFold(n_splits=4))
    - class_weights: class weights to account for class imbalance in performance measurements
    
    Output: average of scores over all CV-runs.
    """
    # define empty matrix to store performances (n CV runs and four performance metrics)
    n_splits_outer = outer_cv.get_n_splits()
    performances = np.zeros((n_splits_outer, 4))
    
    # define outer loop
    loop = 0
    for train_outer, test_outer in outer_cv.split(X, y, groups):
        X_train, X_test = X[train_outer], X[test_outer]
        y_train, y_test = y[train_outer], y[test_outer]
        groups_train, groups_test = groups[train_outer], groups[test_outer]
        
        # define inner loop (in GridSearchCV)
        tuned_model = GridSearchCV(model, cv=inner_cv, param_grid=parameter_grid, scoring=scorer)
        tuned_model.fit(X_train, y_train, groups=groups_train)
        
        # make predictions for test set (outer loop)
        y_pred = tuned_model.predict(X_test)
        
        # evaluate performance (factoring in class imbalance)
        recall_list = list(recall_score(y_test, y_pred, average=None))
        precision_list = list(precision_score(y_test, y_pred, average=None))
        f1_list = list(f1_score(y_test, y_pred, average=None))
        accuracy = accuracy_score(y_test, y_pred)
        recall = sum([a*b for a,b in zip(recall_list, class_weights)])/sum(class_weights)
        precision = sum([a*b for a,b in zip(precision_list, class_weights)])/sum(class_weights)
        f1 = sum([a*b for a,b in zip(f1_list, class_weights)])/sum(class_weights)
        performances[loop,:] = [accuracy, recall, precision, f1]
        
        # next loop
        loop += 1
    
    average_performances = performances.mean(0)
    return average_performances  

def NestedGroupKFoldProba(model, X, y, parameter_grid, groups, class_weights, scorer=make_scorer(accuracy_score), 
                     inner_cv=GroupKFold(n_splits=4), outer_cv=GroupKFold(n_splits=4)):
    """
    Implements a nested version of GroupKFold cross-validation using GridSearchCV to evaluate models 
    that need hyperparameter tuning in settings where different groups exist in the available data.
    
    Dependencies: sklearn.model_selection, sklearn.metrics, numpy
    
    Input:
    - X, y: features and labels (must be NumPy arrays).
    - model, parameter_grid: the model instance and its parameter grid to be optimized.
    - groups: the groups to use in both inner- and outer loop.
    - scorer: the scoring to use in inner loop (default: accuracy).
    - inner_cv, outer_cv: the iterators for both CV-loops (default: GroupKFold(n_splits=4))
    - class_weights: class weights to account for class imbalance in performance measurements
    
    Output: cross-validated predicted class probabilities
    """
    # define empty matrix to store performances (n CV runs and four performance metrics)
    n_classes = len(class_weights)
    probabilities = np.zeros((X.shape[0], n_classes))
    preds = np.zeros(X.shape[0])
    
    # define outer loop
    for train_outer, test_outer in outer_cv.split(X, y, groups):
        X_train, X_test = X[train_outer], X[test_outer]
        y_train, y_test = y[train_outer], y[test_outer]
        groups_train, groups_test = groups[train_outer], groups[test_outer]
        
        # define inner loop (in GridSearchCV)
        tuned_model = GridSearchCV(model, cv=inner_cv, param_grid=parameter_grid, scoring=scorer)
        tuned_model.fit(X_train, y_train, groups=groups_train)
        
        # make predictions for test set (outer loop)
        y_probs = tuned_model.predict_proba(X_test)
        y_preds = tuned_model.predict(X_test)
        
        for i, index in enumerate(test_outer):
            probabilities[index,:] = y_probs[i,:]
            preds[index] = y_preds[i]
    
    return probabilities, preds



# PRECISION AND RECALL IN MULTICLASS SETTING
# --------------------------------------------------
def multiclass_precision_recall(y_true, probas_pred, classes, make_plot=False):
    """
    A multiclass extension of the 'precision_recall_curve' function in Scikit-learn. Also includes a plotting
    functionality and outputs the PR AUC score.
    
    Dependencies: matplotlib.pyplot, sklearn.preprocessing, sklearn.metrics
    
    Input:
    - y-true: true labels
    - probas_pred: multiclass probabilities
    - classes: unique classes needed for binarization (not class weights or counts, see label_binarize in sklearn)
    - make_plot: boolean (default: False)
    
    Output: if plot=False (default): P, R and AUC for every class (dicts); else: a figure.
    """
    # define needed variables
    y_true_binary = label_binarize(y_true, classes)
    n_classes = len(classes)
    precision = dict()
    recall = dict()
    auc_scores = dict()
    
    # calculate P, R and AUC per class
    for i in range(n_classes):
        precision[i], recall[i], _ = precision_recall_curve(y_true_binary[:, i], probas_pred[:, i])
        auc_scores[i] = auc(recall[i], precision[i])
    
    # make plot if asked and return results
    if make_plot:
        fig, ax = plt.subplots(figsize=(10,7))
        for i in range(n_classes):
            ax.plot(recall[i], precision[i], lw=2)
            
        ax.set_xlabel('Recall', size=12)
        ax.set_ylabel('Precision', size=12)
        ax.legend(classes)
        fig.tight_layout()
        return fig
    
    else:
        return precision, recall, auc_scores
    

def multiclass_average_PR(y_true, probas_pred, classes, class_weights):
    """
    Computes a weighted macro avegerage for precision and recall over different thresholds.
    
    Dependencies: numpy
    
    Input:
    - y_true
    - probas_pred: probabilities outputted by classifier
    - classes: unique list of classes (not class weights/counts)
    - class_weights: instance counts per class
    """
    # ignore undefined metric warning
    warnings.filterwarnings(action='ignore', category=UndefinedMetricWarning)

    # binarize & compute pred
    y_true_binary = label_binarize(y_true, classes)
    n_classes = len(classes)
    thresholds = np.arange(0, 1, 0.005)
    P_classes = np.zeros((n_classes, len(thresholds)))
    R_classes = np.zeros((n_classes, len(thresholds)))

    # loop over classes
    for i in range(n_classes):
        y_test = y_true_binary[:,i]
        previous_P = 0
        for j, thres in enumerate(thresholds):
            # select class column and threshold probabilities
            y_pred = (probas_pred[:,i] > thres)*1
            
            # compute P, R
            recall_i = recall_score(y_test, y_pred)
            precision_i = precision_score(y_test, y_pred)
            if (precision_i == 0) & (previous_P == 1):
                precision_i = 1
            previous_P = precision_i
            P_classes[i,j] = precision_i
            R_classes[i,j] = recall_i
        
        # average over classes
        P_average = np.average(P_classes, axis=0, weights=class_weights)
        R_average = np.average(R_classes, axis=0, weights=class_weights)

    return P_average, R_average 


# LOCAL BLAST
# --------------------------------------------------
def blast_local(train_indices, test_indices, database):
    """
    This function performs a local BLASTp for given test_indices against a database consisting of
    all train_inidices of some given database that has a column named 'protein_sequence'. Top hits
    are converted to the corresponding numerical class and finally outputted 
    
    Dependencies: NcbiblastpCommandline (Bio.Blast.Applications), NCBIXML (Bio.Blast)
    """
    # define dictionary of classes
    label_dict_species = {'Staphylococcus aureus': 0, 'Klebsiella pneumoniae': 1, 'Acinetobacter baumannii': 2,
                          'Pseudomonas aeruginosa': 3, 'Escherichia coli': 4, 'Salmonella enterica': 5, 
                          'Clostridium difficile': 6}
    
    # construct database against which to blast
    temp_fasta = open('RBPlogodatabase.fasta', 'w')
    for train_index in train_indices:
        sequence = database['protein_sequence'][train_index]
        host = database['host'][train_index]
        temp_fasta.write('>'+host+'\n'+sequence+'\n')
    temp_fasta.close()
    
    # loop over the sequences of the left out group
    prediction_list = []
    for test_index in test_indices:
        # make temp fasta file for test sequence
        test_sequence = database['protein_sequence'][test_index]
        temp_fasta_test = open('testsequence.fasta', 'w')
        temp_fasta_test.write('>'+str(test_index)+'\n'+test_sequence+'\n')
        temp_fasta_test.close()
        
        # perform actual BLAST
        blastp_cline = NcbiblastpCommandline(query='testsequence.fasta', subject='RBPlogodatabase.fasta', evalue=10, 
                                        outfmt=5, out='resultlogo.xml')
        stdout, stderr = blastp_cline()
        
        # save predictions
        # check title of highest scoring match, title is host
        blast_result = open('resultlogo.xml')
        blast_record = NCBIXML.read(blast_result)
        try:
            prediction = blast_record.descriptions[0].title
            prediction = prediction.split(' ')[1:]
            
            # convert to one of keys
            if len(prediction) > 2:
                prediction = ' '.join(prediction[:2])
            else:
                prediction = ' '.join(prediction)
            prediction_list.append(label_dict_species[prediction])
            
        except IndexError:
            prediction = 'None'
            prediction_list.append(9)

        # remove temp fasta test file
        os.remove('testsequence.fasta')
        os.remove('resultlogo.xml')
    
    # remove temp fasta database
    os.remove('RBPlogodatabase.fasta')
    
    return np.array(prediction_list)


# BLAST PREDICTOR & PLOT
# --------------------------------------------------
def blast_top_prediction(sequence, email, treshold=0.01, top3=True):
    """
    This function collects the top hit(s) of a BLAST search (based on e-value and sequence identity) that is not
    the sequence itself. Then it accesses NCBI to collect the host name related to the top hit.
    
    Input:
        * sequence: protein sequence as string
        * email: email to access NCBI
        * treshold: e-value treshold to consider (default=0.01)
        * top3: return top3 predictions? (default=True)
    """
    
    import re
    import time
    import urllib
    Entrez.email = email
    
    # blast sequence
    connection = 0
    while connection == 0:
        try:
            result = NCBIWWW.qblast('blastp', 'nr', sequence=sequence)
            blast_record = NCBIXML.read(result)
            e_list = []
            host_lst = []
            iden_list = []
            title_list = []
            
            # collect BLAST output
            for i in range(4): # top3, but first sequence might be the query itself...
                alignment = blast_record.alignments[i]
                description = blast_record.descriptions[i]
                evalue = description.e
                hsp = alignment.hsps[0]
                identity = (hsp.identities/hsp.align_length)*100
                e_list.append(evalue)
                iden_list.append(identity)
                title_list.append(description.title)
            connection = 1
        except urllib.error.URLError as err:
            time.sleep(1)
            print('Having connection issues... trying again.')
            pass
        except IndexError:
            time.sleep(1)
            print('IndexError...trying again.')
            pass
            
    
    # define patterns for re matching
    pattern_gi = 'gi.[0-9]+.'
    pattern_acc = '\|.+?\|'
    
    # collect host name(s) from NCBI record
    for j, e in enumerate(e_list):
        if (e < treshold) & (iden_list[j] < 100):
            title = title_list[j]
            match_gi = re.search(pattern_gi, title)
            
            # match gene identifier (gi) or accession number if gi is unavailable
            if match_gi == None:
                match_acc = re.search(pattern_acc, title)
                gi = title[match_acc.start()+1:match_acc.end()-1]
            else:
                gi = title[match_gi.start()+3:match_gi.end()-1]
                      
            # fetch protein
            error = 0
            try:
                handle = Entrez.efetch(db='protein', id=gi, rettype='gb', retmode='text')
            except urllib.error.HTTPError as err:
                print(gi)
                error = 1
                pass
                #if err.code == 400:
            
            # check host
            if error != 1:
                for record in SeqIO.parse(handle, 'genbank'):
                    host = None
                    if 'host' in record.features[0].qualifiers:
                        host = record.features[0].qualifiers['host'][0]
                    elif 'lab_host' in record.features[0].qualifiers:
                        host = record.features[0].qualifiers['lab_host'][0]
                    elif 'strain' in record.features[0].qualifiers:
                        host = record.features[0].qualifiers['strain'][0]
                    elif 'organism' in record.annotations:
                        # parse relevant info with re (everything up to first space)
                        text_org = record.annotations['organism']
                        pattern_org = '.+? '
                        match_org = re.search(pattern_org, text_org)
                        if match_org != None:
                            host = text_org[match_org.start():match_org.end()-1]
                        else:
                            # if organism contains no spaces (e.g. 'Acinetobacter', we can't match with the above pattern...)
                            pattern_org = '.*'
                            match_org = re.search(pattern_org, text_org)
                            host = text_org[match_org.start():match_org.end()]
            
                # cut host name off to species level
                pattern_species = '.+? .+? '
                match_species = re.search(pattern_species, host)
                if match_species is not None:
                    host = host[match_species.start():match_species.end()-1]
            
                # append host to list
                if host != None:
                    host_lst.append(host)
            
    if len(host_lst) == 0:
        print('treshold too strict to predict host(s)')
        return '-'
    elif top3 == False:
        return host_lst[0]
    elif top3 == True:
        if len(host_lst) == 4:
            host_lst = host_lst[0:3]
            return host_lst
        else:
            return host_lst
    
    
def blast_host_predictor(sequence, email, treshold=1e-4):
    """
    This function blasts a given phage protein sequence and returns the hosts of the resulting alignments that have an
    e-value under the given treshold. Put simply, this functions predicts the host(s) of the related phage based 
    on BLAST.
    
    Input: a phage protein sequence, an email for Entrez.efetch and (optionally) a treshold for e-value (default=1e-50)
    Output: a dictionary of hosts and their percentages of occurrence
    
    Remarks: implement extra cutoff for identities (or positives), 
        see blast_record.hsps.identities or .align_length or .positives
    """
    import re
    import time
    import urllib
    Entrez.email = email
    
    # blast sequence
    result = NCBIWWW.qblast('blastp', 'nr', sequence=sequence)
    blast_record = NCBIXML.read(result)
    
    host_lst = []
    pattern_gi = 'gi.[0-9]+.'
    pattern_acc = '\|.+?\|'
    
    # parse blast record for relevant results
    for description in blast_record.descriptions:
        time.sleep(0.15)
        if description.e <= treshold:
            # regular expression to filter gi
            title = description.title
            match_gi = re.search(pattern_gi, title)
            
            # match gene identifier (gi) or accession number if gi is unavailable
            if match_gi == None:
                match_acc = re.search(pattern_acc, title)
                gi = title[match_acc.start()+1:match_acc.end()-1]
            else:
                gi = title[match_gi.start()+3:match_gi.end()-1]
                      
            # fetch protein
            error = 0
            try:
                handle = Entrez.efetch(db='protein', id=gi, rettype='gb', retmode='text')
            except urllib.error.HTTPError as err:
                print(gi)
                error = 1
                pass
                #if err.code == 400:
            
            # check host
            if error != 1:
                for record in SeqIO.parse(handle, 'genbank'):
                    if 'host' in record.features[0].qualifiers:
                        host = record.features[0].qualifiers['host'][0]
                    elif 'lab_host' in record.features[0].qualifiers:
                        host = record.features[0].qualifiers['lab_host'][0]
                    elif 'strain' in record.features[0].qualifiers:
                        host = record.features[0].qualifiers['strain'][0]
                    elif 'organism' in record.annotations:
                        # parse relevant info with re (everything up to first space)
                        text_org = record.annotations['organism']
                        pattern_org = '.+? '
                        match_org = re.search(pattern_org, text_org)
                        if match_org != None:
                            host = text_org[match_org.start():match_org.end()-1]
                        else:
                            # if organism contains no spaces (e.g. 'Acinetobacter', we can't match with the above pattern...)
                            pattern_org = '.*'
                            match_org = re.search(pattern_org, text_org)
                            host = text_org[match_org.start():match_org.end()]
            
                # cut host name off to species level
                pattern_species = '.+? .+? '
                match_species = re.search(pattern_species, host)
                if match_species is not None:
                    host = host[match_species.start():match_species.end()-1]
            
                # append host to list
                host_lst.append(host)
            
    # count hosts (if any)
    if len(host_lst) > 0:
        # count number of times hosts occur 
        host_dict = {}
        for item in host_lst:
            if item in host_dict:
                host_dict[item] += 1
            else:
                host_dict[item] = 1  
        # divide counts to percentages
        for item in host_dict:
            host_dict[item] /= len(host_lst)
        return host_dict
    
    else:
        print('treshold too strict to predict host(s)')
