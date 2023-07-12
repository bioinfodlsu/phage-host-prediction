"""
=======================================================
This script contains the constants used in the project.
@author    Mark Edward M. Gonzales
=======================================================
"""

class ConstantsUtil(object):    
    # ===========
    # Directories
    # ===========
    DATA = 'data'
    FASTA = 'fasta'
    
    HYPOTHETICAL = 'hypothetical'
    RBP = 'rbp'
    NUCLEOTIDE = 'nucleotide'
    COMPLETE = 'complete'

    GENBANK = 'genbank'
    PROKKA = 'prokka'
    MASTER = 'master'
    
    # =============
    # Preprocessing
    # =============
    PREPROCESSING = 'preprocessing'
    GENUS_TYPO = f'{PREPROCESSING}/genus_typo.txt'
    
    BACTERIA_NOT_GENUS = f'{PREPROCESSING}/bacteria_not_genus.txt'
    EXCLUDED_HOSTS = f'{PREPROCESSING}/excluded_hosts.txt'
    
    NCBI_STANDARD_NOMENCLATURE = f'{PREPROCESSING}/ncbi_standard_nomenclature.txt'
        
    # Regex for candidate genera
    CANDIDATE_REGEX = r'candidat(e|us)'
    # Regex for selecting annotated RBPs
    RBP_REGEX = r'tail?(.?|\s*)(?:spike?|fib(?:er|re))|recept(?:o|e)r(.?|\s*)(?:bind|recogn).*(?:protein)?|(?<!\w)RBP(?!a)'
    # Regex for token delimiters in gene product annotations
    TOKEN_DELIMITER = '[-\|,.\/\s]'
    
    HYPOTHETICAL_KEYWORDS = f'{PREPROCESSING}/hypothetical_keywords.txt'
    RBP_RELATED_NOT_RBP = f'{PREPROCESSING}/rbp_related_not_rbp.txt'
    PUTATIVE_FUNCTIONS = f'{PREPROCESSING}/putative_functions.txt'
    
    # Number of entries for displayed progress
    DISPLAY_PROGRESS = 1000
    # Minimum edit distance to be considered a possible misspelling
    MISSPELLING_THRESHOLD = 2
    # Minimum length for a token to be considered a keyword of interest
    MIN_LEN_KEYWORD = 6
    
    # Bounds for the length of RBPs (excluding those with outlying lengths)
    LOWER_BOUND_RBP_LENGTH = -533.0
    UPPER_BOUND_RBP_LENGTH = 1587.0

    # ===============
    # Temporary Files
    # ===============
    TEMP = 'temp'
    TEMP_PREPROCESSING = f'{TEMP}/{PREPROCESSING}'
    RESULTS = 'results'
    TEMP_RESULTS = f'{TEMP}/{RESULTS}'
    
    INPHARED_WITH_HOSTS = 'inphared.csv'
    
    NO_CDS_ANNOT = 'no_cds_annot.pickle'
    ANNOT_PRODUCTS = 'annot_products.pickle'
    RBP_PRODUCTS = 'rbp_products.pickle'
    HYPOTHETICAL_PRODUCTS = 'hypothetical_proteins.pickle'
    
    RBP_LENGTHS = 'rbp_lengths.pickle'
    
    FOR_EMBED = {HYPOTHETICAL: f'{TEMP}/hypothetical_for_embed',
                 RBP: f'{TEMP}/rbp_for_embed'}

    # ========
    # INPHARED
    # ========
    INPHARED = 'inphared'
    INPHARED_RBP_DATA = f'rbp.csv'
    INPHARED_EXTRA_COLS = 36
    
    TEMP_INPHARED = f'{TEMP}/{INPHARED}'
    TOP_GENUS = f'{TEMP_INPHARED}/top_genus.pickle'
    TOP_GENUS_ACCESSION = f'{TEMP_INPHARED}/top_genus_accession.pickle'
    NOT_TOP_GENUS_ACCESSION = f'{TEMP_INPHARED}/not_top_genus_accession.pickle'
    INPHARED_ACCESSION = f'{TEMP_INPHARED}/inphared_accession.pickle'
    
    # =======================
    # Protein Language Models
    # =======================
    EMBEDDINGS = 'embeddings'
    PLM = {'PROTTRANSBERT': f'{EMBEDDINGS}/prottransbert', 
           'PROTXLNET': f'{EMBEDDINGS}/protxlnet',
           'PROTTRANSALBERT': f'{EMBEDDINGS}/prottransalbert',
           'PROTT5': f'{EMBEDDINGS}/prott5', 
           'ESM': f'{EMBEDDINGS}/esm', 
           'ESM1B': f'{EMBEDDINGS}/esm1b', 
           'SEQVEC': f'{EMBEDDINGS}/seqvec'}
    
    EMBEDDINGS_CSV = 'rbp_embeddings'
    PLM_EMBEDDINGS_CSV = {'PROTTRANSBERT': f'{EMBEDDINGS_CSV}_prottransbert.csv', 
                          'PROTXLNET': f'{EMBEDDINGS_CSV}_protxlnet.csv',
                          'PROTTRANSALBERT': f'{EMBEDDINGS_CSV}_prottransalbert.csv',
                          'PROTT5': f'{EMBEDDINGS_CSV}_prott5.csv', 
                          'ESM': f'{EMBEDDINGS_CSV}_esm.csv', 
                          'ESM1B': f'{EMBEDDINGS_CSV}_esm1b.csv', 
                          'SEQVEC': f'{EMBEDDINGS_CSV}_seqvec.csv',
                          'BOECKAERTS': f'{EMBEDDINGS_CSV}_boeckaerts.csv'}
    
    PLM_RESULTS = {'PROTTRANSBERT': f'{TEMP_RESULTS}/prottransbert.pickle', 
                   'PROTXLNET': f'{TEMP_RESULTS}/protxlnet.pickle',
                   'PROTTRANSALBERT': f'{TEMP_RESULTS}/prottransalbert.pickle',
                   'PROTT5': f'{TEMP_RESULTS}/prott5.pickle', 
                   'ESM': f'{TEMP_RESULTS}/esm.pickle', 
                   'ESM1B': f'{TEMP_RESULTS}/esm1b.pickle', 
                   'SEQVEC': f'{TEMP_RESULTS}/seqvec.pickle',
                   'BOECKAERTS': f'{TEMP_RESULTS}/boeckaerts.pickle'}
    
    FEATURE_IMPORTANCE = f'{TEMP}/feature_importance.pickle'
    
    COMPLETE_EMBEDDINGS = f"{INPHARED}/{PLM['PROTTRANSBERT']}/{COMPLETE}/{MASTER}"
    
    TRAINED_MODEL = 'models'
    PLM_TRAINED_MODEL = {'PROTTRANSBERT': f'{TRAINED_MODEL}/prottransbert.joblib', 
                         'PROTXLNET': f'{TRAINED_MODEL}/protxlnet.joblib',
                         'PROTTRANSALBERT': f'{TRAINED_MODEL}/prottransalbert.joblib',
                         'PROTT5': f'{TRAINED_MODEL}/prott5.joblib', 
                         'ESM': f'{TRAINED_MODEL}/esm.joblib', 
                         'ESM1B': f'{TRAINED_MODEL}/esm1b.joblib', 
                         'SEQVEC': f'{TRAINED_MODEL}/seqvec.joblib',
                         'BOECKAERTS': f'{TRAINED_MODEL}/boeckaerts.joblib'}

    # ==============
    # RBP Prediction
    # ==============
    XGB_RBP_PREDICTION = 'rbp_prediction/RBPdetect_xgb_model.json'
    
    def __init__(self, date = ''):
        """
        Constructor
        
        Parameters:
        - date: Download date of the dataset
        """
        self.DATE = date
        self.INPHARED_GENOME = f'/datasets/{self.INPHARED}/{self.INPHARED}/GenomesDB'
        self.INPHARED_TSV = f'{self.INPHARED}/{self.DATE}_data_excluding_refseq.tsv'
        self.INPHARED_GB = f'{self.INPHARED}/{self.DATE}_phages_downloaded_from_genbank.gb'