"""
======================================================================
This script contains the utility functions for sequence preprocessing.
@author    Mark Edward M. Gonzales
======================================================================
"""

import regex as re
import os
from collections import defaultdict

import numpy as np
import nltk
from Bio import SeqIO

class SequenceParsingUtil(object):
    def __init__(self, display_progress = None, misspelling_threshold = None, min_len_keyword = None):
        """
        Constructor
        
        Parameters:
        - display_progress: True if the number of entries processed is to be displayed periodically; False, otherwise
        - misspelling_threshold: Minimum edit distance to be considered a possible misspelling
        - min_len_keyword: Minimum length for a token to be considered a keyword of interest
        """
        self.display_progress = display_progress
        self.misspelling_threshold = misspelling_threshold
        self.min_len_keyword = min_len_keyword
        self.token_delimiter = None
        
        self.inphared_gb = None
        self.inphared = None
        
        self.phages_unspec_host = None
        self.unfiltered_hosts = None
        self.genus_typo = None
        self.unfiltered_suspected_genera = None
        self.excluded_hosts = None
        self.valid_hosts = None
       
        self.hypothetical_keywords = set()
        self.rbp_related_not_rbp = set()
        self.putative_functions = set()
        
        self.no_cds_annot = None
        self.annot_products = None
        self.rbp_products = None
        self.hypothetical_proteins = None
        
    def __display_progress(self, ctr):
        """
        Periodically displays the number of processed entries
        
        Parameters:
        - ctr: Number of processed entries
        """
        if ctr % self.display_progress == 0:
            print("Processed", ctr, "records")
            
    # =======
    # Setters
    # =======
        
    def set_inphared_gb(self, inphared_gb):
        """
        Sets the consolidated GenBank entries of the entries fetched via INPHARED
        
        Parameters:
        - inphared_gb: File path of the consolidated GenBank entries of the entries fetched via INPHARED
        """
        self.inphared_gb = inphared_gb
        
    def set_inphared(self, inphared):
        """
        Sets the dataset (phage-host table, along with other information)
        
        Parameters:
        - inphared: File path of the dataset (phage-host table, along with other information)
        """
        self.inphared = inphared
           
    def set_no_cds_annot(self, no_cds_annot):
        """
        Sets the phage entries without coding sequence information
        
        Parameters:
        - no_cds_annot: Set of phage entries without coding sequence information
        """
        self.no_cds_annot = no_cds_annot
        
    def set_annot_products(self, annot_products):
        """
        Sets the gene product annotations
        
        Parameters:
        - annot_products: Set of gene product annotations
        """
        self.annot_products = annot_products
        
    def set_rbp_products(self, rbp_products):
        """
        Sets the RBP annotations
        
        Parameters:
        - rbp_products: Set of RBP annotations
        """
        self.rbp_products = rbp_products
        
    def set_hypothetical_proteins(self, hypothetical_proteins):
        """
        Sets the hypothetical protein annotations
        
        Parameters:
        - hypothetical_proteins: Set of hypothetical protein annotations
        """
        self.hypothetical_proteins = hypothetical_proteins
        
    def set_token_delimiter(self, token_delimiter):
        """
        Sets the possible token delimiters for gene product annotations
        
        Parameters:
        - token_delimiters: Regular expression for the possible token delimiters for gene product annotations
        """
        self.token_delimiter = token_delimiter

    # ==================
    # Data Preprocessing
    # ==================
        
    def get_phages_unspec_host(self, inphared_unspec_host):
        """
        Retrieves the isolation host information for phages where the returned host information of INPHARED is "unspecified"
        
        Parameters:
        - inphared_unspec_host: DataFrame containing the phage entries with unspecified host information
        
        Returns:
        - Dictionary where each key is a phage with unspecified host information and the value is the list of isolation
          hosts retrieved from GenBank
        """
        records = SeqIO.parse(self.inphared_gb, 'gb')

        accession_np = inphared_unspec_host['Accession'].to_numpy()
        phages_unspec_host = defaultdict(list)

        ctr = 0
        iter_flag = True
        while iter_flag:
            try:
                record = next(records)

                try:            
                    idx = np.where(accession_np == record.name)[0][0]

                    if 'host' in record.features[0].qualifiers:
                        phages_unspec_host[record.name].append(record.features[0].qualifiers['host'])
                    else:
                        phages_unspec_host[record.name].append([])

                    if 'lab_host' in record.features[0].qualifiers:
                        phages_unspec_host[record.name].append(record.features[0].qualifiers['lab_host'])
                    else:
                        phages_unspec_host[record.name].append([])

                except IndexError:
                    pass

                ctr += 1
                self.__display_progress(ctr)

            except StopIteration:
                iter_flag = False
        
        self.phages_unspec_host = phages_unspec_host
        
        return phages_unspec_host
    
    def get_unfiltered_hosts(self):
        """
        Retrieves the isolation hosts of phages where the returned host information of INPHARED is "unspecified"
        
        Returns:
        - Set of isolation hosts of phages where the returned host information of INPHARED is "unspecified"
        """
        unfiltered_hosts = set()
        for key, value in self.phages_unspec_host.items():
            for host in value[0]:
                unfiltered_hosts.add(host)
            for host in value[1]:
                unfiltered_hosts.add(host)
                
        self.unfiltered_hosts = unfiltered_hosts

        return unfiltered_hosts
    
    def get_genus_typo(self, filename):
        """
        Retrieves the genera with typographical errors
        
        Parameters:
        - filename: File path of the text file recording the genera with typographical errors
        
        Returns:
        - Dictionary where each key is the misspelled genus and the value is the correct spelling
        """
        genus_typo = {}
        with open(filename) as typo_file:
            for entry in typo_file:
                typo, correct = entry.split('\t')
                typo = typo.rstrip('\n')
                correct = correct.rstrip('\n')
                genus_typo[typo] = correct
                
        self.genus_typo = genus_typo
                
        return genus_typo
    
    def __get_suspected_genus(self, candidate_regex, host):
        """
        Returns the genus of the isolation host of a phage where the returned host information of INPHARED is "unspecified"
        
        Parameters:
        - candidate_regex: Regex for matching candidate genera
        
        Returns:
        - Genus of the isolation host of a phage where the returned host information of INPHARED is "unspecified"
        """
        host_tokens = host.split(' ')

        if re.search(candidate_regex, host_tokens[0], re.IGNORECASE):
            # Handle the case where name starts with "Candidate division"
            if host_tokens[1] == 'division':
                genus = host_tokens[2].lower()
            else:
                genus = host_tokens[1].lower()
        else:
            genus = host_tokens[0].lower()

        if genus in self.genus_typo:
            genus = self.genus_typo[genus]

        return genus
    
    def get_unfiltered_suspected_genera(self, candidate_regex):
        """
        Returns the genera of the isolation hosts of phages where the returned host information of INPHARED is "unspecified"
        
        Parameters:
        - candidate_regex: Regex for matching candidate genera
        
        Returns:
        - Set of genera of the isolation hosts of phages where the returned host information of INPHARED is "unspecified"
        """
        unfiltered_suspected_genera = set()

        for unfiltered_host in self.unfiltered_hosts:
            unfiltered_suspected_genera.add(self.__get_suspected_genus(candidate_regex, unfiltered_host))
            
        self.unfiltered_suspected_genera = unfiltered_suspected_genera

        return unfiltered_suspected_genera
    
    def get_excluded_hosts(self, *exclusion_files):
        """
        Retrieves the isolation hosts that do not pertain to bacterial hosts or those that are less specific than genus level
        
        Parameters:
        - exclusion files: File path(s) of the text file(s) recording the isolation hosts in GenBank that do not pertain
          to bacterial hosts or those that are less specific than genus level
        
        Returns:
        - Set of isolation hosts that do not pertain to bacterial hosts or those that are less specific than genus level
        """
        excluded_hosts = set()
        
        for exclusion_file in exclusion_files:
            with open(exclusion_file, 'r') as file:
                for host in file:
                    excluded_hosts.add(host.rstrip('\n'))
        
        self.excluded_hosts = excluded_hosts
        
        return excluded_hosts
    
    def get_valid_hosts(self):
        """
        Returns the isolation hosts that pertain to bacterial hosts at the genus level
        
        Returns:
        - Set of isolation hosts that pertain to bacterial hosts at the genus level
        """
        self.valid_hosts = self.unfiltered_suspected_genera - self.excluded_hosts
        return self.valid_hosts
    
    def is_possible_misspelling(self, typo, keyword):
        """
        Returns True if a word is a possible misspelling; False, otherwise. This function uses the minimum edit distance
        as a heuristic for identifying possible misspellings
        
        Parameters:
        - typo: Potentially misspelled word
        - keyword: Correct spelling
        
        Returns:
        - True if a word is a possible misspelling; False, otherwise. This function uses the minimum edit distance
          as a heuristic for identifying possible misspellings
        """
        return (len(keyword) >= self.min_len_keyword and 
                nltk.edit_distance(typo, keyword, transpositions = True) <= self.misspelling_threshold)
    
    def update_host_column(self, candidate_regex, inphared_unspec_host, inphared_augmented):
        """
        Updates the host column of the phage-host table with the isolation host information from GenBank
        
        Parameters:
        - candidate_regex: Regex for matching candidate genera
        - inphared_unspec_host: DataFrame containing the phage entries with unspecified host information
        - inphared_augmented: Original phage-host table to be updated
        """
        records = SeqIO.parse(self.inphared_gb, 'gb')

        accession_unspec_np = inphared_unspec_host['Accession'].to_numpy()
        accession_augmented_np = inphared_augmented['Accession'].to_numpy()

        ctr = 0
        iter_flag = True
        while iter_flag:
            try:
                record = next(records)

                try:
                    unspec_idx = np.where(accession_unspec_np == record.name)[0][0]

                    hosts = set()

                    if 'host' in record.features[0].qualifiers:
                        for host in record.features[0].qualifiers['host']:
                            host_genus = self.__get_suspected_genus(candidate_regex, host)
                            if host_genus in self.valid_hosts:
                                hosts.add(host_genus)

                    if 'lab_host' in record.features[0].qualifiers:
                        for host in record.features[0].qualifiers['lab_host']:
                            host_genus = self.__get_suspected_genus(candidate_regex, host)
                            if host_genus in self.valid_hosts:
                                hosts.add(host_genus)

                    if hosts:
                        hosts_str = ' | '.join(hosts)
                        augmented_idx = np.where(accession_augmented_np == record.name)[0][0]
                        inphared_augmented.at[augmented_idx, 'Host'] = hosts_str

                        print("Index:", augmented_idx, "| Name:", record.name, "| Host:", hosts_str)

                except IndexError:
                    pass

                # Display progress
                ctr += 1
                self.__display_progress(ctr)

            except StopIteration:
                iter_flag = False
    
    def get_ncbi_standard_nomenclature(self, filename):
        """
        Retrieves the standard nomenclature following NCBI Taxonomy
        
        Parameters:
        - filename: File path of the text file containing the equivalent nomenclature of a genus following NCBI Taxonomy
        """
        nomenclature = {}
        with open(filename) as nomenclature_file:
            for entry in nomenclature_file:
                nonstandard, standard = entry.split('\t')
                nonstandard = nonstandard.rstrip('\n')
                standard = standard.rstrip('\n')
                nomenclature[nonstandard] = standard
                                
        return nomenclature
    
    # ==================
    # CDS Identification
    # ==================
    
    def get_no_cds_annot(self):
        """
        Returns the set of phage entries without coding sequence information
        
        Returns:
        - Set of phage entries without coding sequence information
        """
        records = SeqIO.parse(self.inphared_gb, 'gb')
        
        accession_np = self.inphared['Accession'].to_numpy()

        ctr = 0
        iter_flag = True

        no_cds_annot = set()
        while iter_flag:
            try:
                record = next(records)

                try:
                    idx = np.where(self.inphared == record.name)[0][0]

                    has_cds_annot = False
                    for feature in record.features:
                        if feature.type == 'CDS':
                            has_cds_annot = True
                            break

                    if not has_cds_annot:
                        no_cds_annot.add((idx, record.name))

                except IndexError:
                    pass

                # Display progress
                ctr += 1
                self.__display_progress(ctr)
                
            except StopIteration:
                iter_flag = False
                
        self.no_cds_annot = no_cds_annot
        
        return no_cds_annot
    
    # ====================================
    # Process entries with CDS information
    # ====================================
    
    def construct_keyword_list(self, hypothetical_file, rbp_related_file, putative_functions_file):
        """
        Constructs a list of keywords associated with hypothetical proteins, proteins with putative functions, and RBP-related
        proteins that are not RBPs themselves
        
        Parameters:
        - hypothetical_file: File path of the text file containing keywords associated with hypothetical proteins
        - rbp_related_file: File path of the text file containing keywords associated with RBP-related proteins 
                            that are not RBPs themselves
        - putative_functions_file: File path of the text file containing keywords associated with proteins with
                                   putative functions
        """
        with open(hypothetical_file) as file:
            for keyword in file:
                self.hypothetical_keywords.add(keyword.rstrip('\n'))

        with open(rbp_related_file) as file:
            for keyword in file:
                self.rbp_related_not_rbp.add(keyword.rstrip('\n'))

        with open(putative_functions_file) as file:
            for keyword in file:
                self.putative_functions.add(keyword.rstrip('\n'))

        self.putative_functions = self.putative_functions.union(self.rbp_related_not_rbp)
    
    def get_annot_products(self):
        """
        Retrieves the annotations for the gene products in GenBank
        
        Returns:
        - Set of annotations for the gene products in GenBank
        """
        records = SeqIO.parse(self.inphared_gb, 'gb')

        accession_np = self.inphared['Accession'].to_numpy()

        ctr = 0
        iter_flag = True

        annot_products = set()
        while iter_flag:
            try:
                record = next(records)

                try:
                    idx = np.where(self.inphared == record.name)[0][0]

                    has_cds_annot = False
                    for feature in record.features:
                        if feature.type == 'CDS':
                            try:
                                annot_products.add(feature.qualifiers['product'][0])
                            except KeyError:
                                pass

                except IndexError:
                    pass

                # Display progress
                ctr += 1
                self.__display_progress(ctr)

            except StopIteration:
                iter_flag = False
                
        self.annot_products = annot_products
        
        return annot_products
    
    def __is_rbp_related_but_not_rbp(self, product):
        """
        Checks if a given annotation pertains to an RBP-related protein that is not an RBP itself
        
        Parameters:
        - product: Annotation for a gene product
        
        Returns:
        - True if a given annotation pertains to an RBP-related protein that is not an RBP itself; False, otherwise
        """
        for keyword in self.rbp_related_not_rbp:
            if keyword in product:
                return True
            else:
                for token in re.split(self.token_delimiter, product):
                    # 'unamed' and 'named' are within the minimum edit distance threshold for misspellings if token is 'rapid'
                    if token != 'rapid' and self.is_possible_misspelling(token, keyword):
                        return True

        return False
    
    def __has_putative_function(self, product):
        """
        Checks if a given annotation pertains to a protein with a putative function
        
        Parameters:
        - product: Annotation for a gene product
        
        Returns:
        - True if a given annotation pertains to a protein with a putative function; False, otherwise
        """
        for keyword in self.putative_functions:
            if keyword in product:
                return True
            else:
                for token in re.split(self.token_delimiter, product):
                    # 'unamed' and 'named' are within the minimum edit distance threshold for misspellings if token is 'rapid'
                    if token != 'rapid' and self.is_possible_misspelling(token, keyword):
                        return True
                    
                    # Enzymes usually end with '-ase'
                    if 'ase' in token:
                        return True                        

        return False
    
    def get_rbp_hypothetical_proteins(self, rbp_regex):
        """
        Construct a list of annotations for RBPs and hypothetical proteins from GenBank annotations
        
        Parameters:
        - rbp_regex: Regex for selecting annotated RBPs
        
        Returns:
        - List of annotations for RBPs from GenBank annotations
        - List of annotations for hypothetical proteins from GenBank annotations
        """
        rbp_products = set()
        hypothetical_proteins = set()

        ctr = 0
        for annot_product in self.annot_products:
            annot_product_lower = annot_product.lower()

            if not self.__is_rbp_related_but_not_rbp(annot_product_lower) and re.search(rbp_regex, annot_product_lower, re.IGNORECASE):
                rbp_products.add(annot_product_lower)

            else:
                hypothetical_keyword_found = False

                if not self.__has_putative_function(annot_product_lower):
                    for keyword in self.hypothetical_keywords:
                        if not hypothetical_keyword_found:
                            if keyword in annot_product_lower:
                                hypothetical_proteins.add(annot_product_lower)
                                hypothetical_keyword_found = True
                            else:
                                tokens = re.split(self.token_delimiter, annot_product_lower)
                                for token in tokens:
                                    # 'unamed' and 'named' are within the minimum edit distance threshold for misspellings
                                    if token != 'named' and self.is_possible_misspelling(token, keyword):
                                        hypothetical_proteins.add(annot_product_lower)
                                        hypothetical_keyword_found = True
                                        break
                        else:
                            break

            ctr += 1
            self.__display_progress(ctr)
                
        self.rbp_products = rbp_products
        self.hypothetical_proteins = hypothetical_proteins
        
        return rbp_products, hypothetical_proteins
    
    # ===================================
    # Analyze distribution of RBP lengths
    # ===================================
    
    def __has_unknown_amino_acid(self, sequence):
        """
        Checks if a given protein sequence has an undetermined or unrecognized amino acid
        
        Parameters:
        - sequence: Protein sequence to be checked
        
        Returns:
        - True if the protein sequence has an undetermined or unrecognized amino acid; False, otherwise
        """
        unknown_aa_expr = r'[^ACDEFGHIKLMNPQRSTVWY]'
        if re.search(unknown_aa_expr, sequence):
            return True

        return False
 
    def __is_acceptable_length(self, sequence, lower_bound, upper_bound):
        """
        Checks if a given protein sequence does not have an outlying length
        
        Parameters:
        - sequence: Protein sequence to be checked
        - lower_bound: Lower bound of the lengths of RBPs (excluding outliers)
        - upper_bound: Upper bound of the lengths of RBPs (excluding outliers)
        
        Returns:
        - True if the protein sequence does not have an outlying length; False, otherwise
        """
        return lower_bound <= len(sequence) and len(sequence) <= upper_bound
    
    def __is_rbp(self, product, sequence, lower_bound = -1, upper_bound = -1):
        """
        Checks if a given gene product is an RBP
        
        Parameters:
        - product: Annotation for the gene product
        - sequence: Protein sequence of the gene product
        - lower_bound: Lower bound of the lengths of RBPs (excluding outliers)
        - upper_bound: Upper bound of the lengths of RBPs (excluding outliers)
        
        Returns:
        - True if the gene product is an RBP; False, otherwise
        """
        if lower_bound == -1 and upper_bound == -1:
            return (product in self.rbp_products
                    and not self.__has_unknown_amino_acid(sequence))
        
        return (product in self.rbp_products
                and not self.__has_unknown_amino_acid(sequence)
                and self.__is_acceptable_length(sequence, lower_bound, upper_bound))
    
    def generate_rbp_len_distribution(self):
        """
        Returns a dictionary containing the counts of the lengths (number of amino acids) of the RBPs based on GenBank annotation
        
        Returns:
        - Dictionary where each key is an RBP length and the value is its count (i.e., the number of RBPs with that length)
        """
        records = SeqIO.parse(self.inphared_gb, 'gb')

        accession_np = self.inphared['Accession'].to_numpy()

        ctr = 0
        iter_flag = True

        len_distribution = defaultdict(lambda: 0)
        while iter_flag:
            try:
                record = next(records)

                try:
                    idx = np.where(self.inphared == record.name)[0][0]
                    name = record.name

                    rbp_fasta_str = ''
                    hypothetical_fasta_str = ''

                    for feature in record.features:
                        if feature.type == 'CDS':
                            try:
                                product = feature.qualifiers['product'][0].lower()
                                translation = feature.qualifiers['translation'][0]

                                if self.__is_rbp(product, translation):
                                    len_distribution[len(translation)] += 1

                            except KeyError:
                                pass

                except IndexError:
                    pass

                # Display progress
                ctr += 1
                self.__display_progress(ctr)

            except StopIteration:
                iter_flag = False
                
        return len_distribution
    
    def generate_rbp_len_distribution_prokka(self, len_distribution, complete_genome_dir):
        """
        Updates the dictionary containing the counts of the lengths of the RBPs with the lengths of the RBPs annotated 
        using Prokka
        
        Parameters:
        - len_distribution: Dictionary containing the counts of the lengths of the RBPs based on GenBank annotation
        - complete_genome_dir: File path of the directory with the genomes of all the phage entries retrieved via INPHARED
        """
        ctr = 0
        for phage in self.no_cds_annot:
            with open(os.path.dirname(os.getcwd()) + f'/{complete_genome_dir}/{phage[1]}/{phage[1]}.faa') as handle:                
                    
                for record in SeqIO.parse(handle, 'fasta'):
                    product = record.description.split(' ', 1)[1]
                    translation = str(record.seq)
                                        
                    if self.__is_rbp(product, translation):
                        len_distribution[len(translation)] += 1

            # Display progress
            ctr += 1
            self.__display_progress(ctr)
    
    # ===============================================
    # Generate FASTA for entries with CDS information
    # ===============================================    
    def __is_hypothetical_protein(self, product, sequence, lower_bound, upper_bound):
        """
        Checks if a given gene product is a hypothetical protein with length within the bounds for RBP lengths
        
        Parameters:
        - product: Annotation for the gene product
        - sequence: Protein sequence of the gene product
        - lower_bound: Lower bound of the lengths of RBPs (excluding outliers)
        - upper_bound: Upper bound of the lengths of RBPs (excluding outliers)
        
        Returns:
        - True if the gene product is a hypothetical protein with length less within the bounds for RBP lengths
        """
        return (product in self.hypothetical_proteins 
                and not self.__has_unknown_amino_acid(sequence) 
                and self.__is_acceptable_length(sequence, lower_bound, upper_bound))
    
    def generate_rbp_hypothetical_fasta(self, rbp_genbank_dir, hypothetical_genbank_dir, lower_bound, upper_bound):
        """
        Generates the FASTA files containing the proteomes of the phages with coding sequence information in GenBank
        
        Parameters:
        - rbp_genbank_dir: File path of the directory with the RBP protein sequences of the phages with coding sequence
                           information in GenBank
        - hypothetical_genbank_dir: File path of the directory with the hypothetical protein sequences of the phages 
                                    with coding sequence information in GenBank
        - lower_bound: Lower bound of the lengths of RBPs (excluding outliers)
        - upper_bound: Upper bound of the lengths of RBPs (excluding outliers)
        """
        records = SeqIO.parse(self.inphared_gb, 'gb')

        accession_np = self.inphared['Accession'].to_numpy()

        ctr = 0
        iter_flag = True

        while iter_flag:
            try:
                record = next(records)

                try:
                    idx = np.where(self.inphared == record.name)[0][0]
                    name = record.name

                    rbp_fasta_str = ''
                    hypothetical_fasta_str = ''

                    for feature in record.features:
                        if feature.type == 'CDS':
                            try:
                                product = feature.qualifiers['product'][0].lower()
                                translation = feature.qualifiers['translation'][0]

                                if self.__is_rbp(product, translation, lower_bound, upper_bound):
                                    protein_id = feature.qualifiers['protein_id'][0]    
                                    rbp_fasta_str += f'>{protein_id} {product} \n{translation}\n'

                                elif self.__is_hypothetical_protein(product, translation, lower_bound, upper_bound):
                                    protein_id = feature.qualifiers['protein_id'][0]    
                                    hypothetical_fasta_str += f'>{protein_id} {product} \n{translation}\n'

                            except KeyError:
                                pass

                    if len(rbp_fasta_str) > 0:
                        file_name = f'{name}-rbp.fasta'

                        with open(os.path.join(rbp_genbank_dir, file_name), 'w') as rbp_fasta_file:
                            rbp_fasta_file.write(rbp_fasta_str)

                    if len(hypothetical_fasta_str) > 0:
                        file_name = f'{name}-hypothetical.fasta'

                        with open(os.path.join(hypothetical_genbank_dir, file_name), 'w') as hypothetical_fasta_file:
                            hypothetical_fasta_file.write(hypothetical_fasta_str)

                except IndexError:
                    pass

                # Display progress
                ctr += 1
                self.__display_progress(ctr)

            except StopIteration:
                iter_flag = False
                                 
    def check_fasta_embeddings_per_phage(self, suffix, fasta_dir, embed_dir):
        """
        Returns the difference between the phages in the embeddings and in the FASTA directories.
        There should be a one-to-one correspondence between the phages in the embeddings and in the FASTA directories
        
        Parameters:
        - suffix: 'hypothetical' for hypothetical proteins or 'rbp' for RBPs
        - fasta_dir: File path of the directory containing the FASTA files with the RBP sequences
        - embed_dir: File path of the directory containing the embeddings of the RBP sequences
        
        Returns:
        - Set of phages in the embeddings directory but not in the FASTA directory
        - Set of phages in the FASTA directory but not in the embeddings directory
        """
        fasta = set()
        for file in os.listdir(fasta_dir):
            fasta.add(file[:-(1 + len(suffix) + len('.fasta'))])

        embeddings = set()
        for file in os.listdir(embed_dir):
            embeddings.add(file[:-(1 + len(suffix) + len('-embeddings.csv'))])

        return embeddings - fasta, fasta - embeddings
    
    def check_fasta_embeddings_per_protein(self, suffix, fasta_dir, embed_dir):
        """
        Prints the difference between the protein IDs in the embeddings and in the FASTA directories.
        There should be a one-to-one correspondence between the protein IDs in the embeddings and in the FASTA directories
        
        Parameters:
        - suffix: 'hypothetical' for hypothetical proteins or 'rbp' for RBPs
        - fasta_dir: File path of the directory containing the FASTA files with the RBP sequences
        - embed_dir: File path of the directory containing the embeddings of the RBP sequences
        
        Returns:
        - List of phages where there are differences in the protein IDs in the embeddings and in the FASTA directories
        """
        erroneous = []
        
        ctr = 0
        for file in os.listdir(fasta_dir):
            with open(f'{fasta_dir}/{file}', 'r') as fasta:                
                proteins_fasta, num_proteins_fasta = self.__get_proteins_fasta(fasta)
            
            phage_name = file[:-(1 + len(suffix) + len('.fasta'))]
            with open(f'{embed_dir}/{phage_name}-{suffix}-embeddings.csv', 'r') as embeddings:
                proteins_embed, num_proteins_embed = self.get_proteins_embeddings(embeddings)
            
            if proteins_fasta != proteins_embed or num_proteins_fasta != num_proteins_embed:
                print(phage_name, '|', proteins_fasta - proteins_embed, '|', proteins_embed - proteins_fasta, '|',
                      num_proteins_fasta, '|', num_proteins_embed)
                erroneous.append(phage_name)
                
            # Display progress
            ctr += 1
            self.__display_progress(ctr)
                        
        return erroneous
    
    def __get_proteins_fasta(self, fasta):
        """
        Returns the protein IDs and the number of proteins in a file containing the proteome of a phage
        
        Parameters:
        - fasta: Contents of the file containing the proteome of a phage
        
        Returns:
        - Set of protein IDs in the file containing the proteome of a phage
        - Number of proteins in the file containing the proteome of a phage
        """
        proteins = set()
        num_lines = 0
        
        is_comment = True
        for line in fasta:
            if is_comment:
                protein_name = line.split(' ')[0][1:]
                proteins.add(protein_name)
                
                num_lines += 1
                
            is_comment = not is_comment
            
        return proteins, num_lines
    
    def get_proteins_embeddings(self, embeddings):
        """
        Returns the protein IDs and the number of proteins in a file containing the embeddings of the protein sequences
        of a phage
        
        Parameters:
        - embeddings: Contents of the file containing the embeddings of the protein sequences of a phage
        
        Returns:
        - Set of protein IDs in the file containing the embeddings of the protein sequences of a phage
        - Number of proteins in the file containing the embeddings of the protein sequences of a phage
        """
        proteins = set()
        # Ignore header row
        num_lines = -1
        
        for line in embeddings:
            protein_name = line.split(',')[0]
            proteins.add(protein_name)
            
            num_lines += 1
            
        proteins.remove('ID')
        
        return proteins, num_lines
        
    # =======================================
    # Process entries without CDS information
    # =======================================    
    
    def get_annot_products_prokka(self, complete_genome_dir):
        """
        Retrieves the annotations for the gene products predicted using PHROG
        
        Parameters:
        - complete_genome_dir: File path of the directory with the genomes of all the phage entries retrieved via INPHARED
        
        Returns:
        - Set of annotations for the gene products predicted using PHROG
        """
        annot_products = set()
        phages = set()
        for folder in os.listdir(os.path.dirname(os.getcwd()) + f'/{complete_genome_dir}'):
            if folder[-1] != 'a':
                phages.add(folder)
        
        ctr = 0
        for phage in phages:
            try:
                with open(os.path.dirname(os.getcwd()) + f'/{complete_genome_dir}/{phage}/{phage}.faa') as handle:
                    for record in SeqIO.parse(handle, 'fasta'):
                        annot_products.add(record.description.split(' ', 1)[1])
            except FileNotFoundError:
                continue
                    
            ctr += 1
            self.__display_progress(ctr)
                                    
        return annot_products
    
    def get_rbp_hypothetical_proteins_prokka(self, rbp_regex, annot_products):
        """
        Construct a list of annotations for RBPs and hypothetical proteins from PHROG annotations 
        
        Parameters:
        - rbp_regex: Regex for selecting annotated RBPs
        - annot_products: Set of annotations for the gene products predicted using PHROG
        
        Returns:
        - List of annotations for RBPs from PHROG annotations
        - List of annotations for hypothetical proteins from PHROG annotations
        """
        rbp_products = set()
        hypothetical_proteins = set()

        ctr = 0
        for annot_product in annot_products:
            annot_product_lower = annot_product.lower()

            if not self.__is_rbp_related_but_not_rbp(annot_product_lower) and re.search(rbp_regex, annot_product_lower, re.IGNORECASE):
                rbp_products.add(annot_product_lower)

            else:
                hypothetical_keyword_found = False

                if not self.__has_putative_function(annot_product_lower):
                    for keyword in self.hypothetical_keywords:
                        if not hypothetical_keyword_found:
                            if keyword in annot_product_lower:
                                hypothetical_proteins.add(annot_product_lower)
                                hypothetical_keyword_found = True
                            else:
                                tokens = re.split(self.token_delimiter, annot_product_lower)
                                for token in tokens:
                                    # 'unamed' and 'named' are within the minimum edit distance threshold for misspellings
                                    if token != 'named' and self.is_possible_misspelling(token, keyword):
                                        hypothetical_proteins.add(annot_product_lower)
                                        hypothetical_keyword_found = True
                                        break
                        else:
                            break

            ctr += 1
            self.__display_progress(ctr)
        
        return rbp_products, hypothetical_proteins
    
    def generate_rbp_hypothetical_fasta_prokka(self, complete_genome_dir, rbp_prokka_dir, hypothetical_prokka_dir, 
                                               lower_bound, upper_bound):
        """
        Generates the FASTA files containing the proteomes of the phages annotated using Prokka
        
        Parameters:
        - rbp_genbank_dir: File path of the directory with the RBP protein sequences of the phages annotated using Prokka
        - hypothetical_genbank_dir: File path of the directory with the hypothetical protein sequences of the phages 
                                    annotated using Prokka
        - lower_bound: Lower bound of the lengths of RBPs (excluding outliers)
        - upper_bound: Upper bound of the lengths of RBPs (excluding outliers)
        """
        ctr = 0
        for phage in self.no_cds_annot:
            with open(os.path.dirname(os.getcwd()) + f'/{complete_genome_dir}/{phage[1]}/{phage[1]}.faa') as handle:                
                rbp_fasta_str = ''
                hypothetical_fasta_str = ''
                    
                for record in SeqIO.parse(handle, 'fasta'):
                    product = record.description.split(' ', 1)[1]
                    translation = str(record.seq)
                                        
                    if self.__is_rbp(product, translation, lower_bound, upper_bound):
                        protein_id = record.name 
                        rbp_fasta_str += f'>{protein_id} {product} \n{translation}\n'

                    elif self.__is_hypothetical_protein(product, translation, lower_bound, upper_bound):
                        protein_id = record.name
                        hypothetical_fasta_str += f'>{protein_id} {product} \n{translation}\n'

                    if len(rbp_fasta_str) > 0:
                        file_name = f'{phage[1]}-rbp.fasta'

                        with open(os.path.join(rbp_prokka_dir, file_name), 'w') as rbp_fasta_file:
                            rbp_fasta_file.write(rbp_fasta_str)
                            
                    if len(hypothetical_fasta_str) > 0:
                        file_name = f'{phage[1]}-hypothetical.fasta'

                        with open(os.path.join(hypothetical_prokka_dir, file_name), 'w') as hypothetical_fasta_file:
                            hypothetical_fasta_file.write(hypothetical_fasta_str)
                            
            # Display progress
            ctr += 1
            self.__display_progress(ctr)
    
    # =============================
    # Generate nucleotide sequences
    # =============================    
                
    def generate_rbp_hypothetical_nucleotide(self, rbp_genbank_dir, hypothetical_genbank_dir, lower_bound, upper_bound):
        """
        Generates the FFN files containing the genomes of the phages with coding sequence information in GenBank
        
        Parameters:
        - rbp_genbank_dir: File path of the directory with the RBP protein sequences of the phages with coding sequence
                           information in GenBank
        - hypothetical_genbank_dir: File path of the directory with the hypothetical protein sequences of the phages with
                                    coding sequence information in GenBank
        - lower_bound: Lower bound of the lengths of RBPs (excluding outliers)
        - upper_bound: Upper bound of the lengths of RBPs (excluding outliers)
        """
        records = SeqIO.parse(self.inphared_gb, 'gb')

        accession_np = self.inphared['Accession'].to_numpy()

        ctr = 0
        iter_flag = True

        while iter_flag:
            try:
                record = next(records)

                try:
                    idx = np.where(self.inphared == record.name)[0][0]
                    name = record.name

                    rbp_fasta_str = ''
                    hypothetical_fasta_str = ''

                    for feature in record.features:
                        if feature.type == 'CDS':
                            try:
                                product = feature.qualifiers['product'][0].lower()
                                translation = feature.qualifiers['translation'][0]

                                if self.__is_rbp(product, translation, lower_bound, upper_bound):
                                    protein_id = feature.qualifiers['protein_id'][0]    
                                    rbp_fasta_str += f'>{protein_id} {product} \n{feature.extract(record.seq)}\n'

                                elif self.__is_hypothetical_protein(product, translation, lower_bound, upper_bound):
                                    protein_id = feature.qualifiers['protein_id'][0]    
                                    hypothetical_fasta_str += f'>{protein_id} {product} \n{feature.extract(record.seq)}\n'

                            except KeyError:
                                pass

                    if len(rbp_fasta_str) > 0:
                        file_name = f'{name}-rbp.ffn'

                        with open(os.path.join(rbp_genbank_dir, file_name), 'w') as rbp_fasta_file:
                            rbp_fasta_file.write(rbp_fasta_str)

                    if len(hypothetical_fasta_str) > 0:
                        file_name = f'{name}-hypothetical.ffn'

                        with open(os.path.join(hypothetical_genbank_dir, file_name), 'w') as hypothetical_fasta_file:
                            hypothetical_fasta_file.write(hypothetical_fasta_str)

                except IndexError:
                    pass

                # Display progress
                ctr += 1
                self.__display_progress(ctr)

            except StopIteration:
                iter_flag = False
                
    
    def get_seq_in_fasta(self, protein_id, fasta):
        """
        Gets the sequence associated with a given protein ID
        
        Parameters:
        - protein_id: Protein ID
        - fasta: File path of the FASTA file with the protein sequences
        
        Returns:
        - Sequence associated with the protein ID
        """
        for record in SeqIO.parse(fasta, 'fasta'):
            record_id = record.description.split(' ', 1)[0]
            
            if record_id == protein_id:
                return str(record.seq)
            
        return None
    
    def generate_rbp_hypothetical_nucleotide_prokka(self, complete_genome_dir, rbp_prokka_dir, hypothetical_prokka_dir,
                                                    lower_bound, upper_bound):
        """
        Generates the FFN files containing the genomes of the phages annotated using Prokka
        
        Parameters:
        - complete_genome_dir: File path of the directory with the genomes of all the phage entries retrieved via INPHARED
        - rbp_prokka_dir: File path of the directory with the RBP protein sequences of the phages annotated using Prokka
        - hypothetical_prokka_dir: File path of the directory with the hypothetical protein sequences of the phages
                                   annotated using Prokka
        - lower_bound: Lower bound of the lengths of RBPs (excluding outliers)
        - upper_bound: Upper bound of the lengths of RBPs (excluding outliers)
        """
        ctr = 0
        for phage in self.no_cds_annot:
            with open(os.path.dirname(os.getcwd()) + f'{complete_genome_dir}/{phage[1]}/{phage[1]}.faa') as handle:                
                rbp_fasta_str = ''
                hypothetical_fasta_str = ''
                    
                for record in SeqIO.parse(handle, 'fasta'):
                    product = record.description.split(' ', 1)[1]
                    translation = str(record.seq)
                                        
                    if self.__is_rbp(product, translation, lower_bound, upper_bound):
                        protein_id = record.name 
                        fasta = os.path.dirname(os.getcwd()) + f'{complete_genome_dir}/{phage[1]}/{phage[1]}.ffn'
                        rbp_fasta_str += f'>{protein_id} {product} \n{self.get_seq_in_fasta(protein_id, fasta)}\n'

                    elif self.__is_hypothetical_protein(product, translation, lower_bound, upper_bound):
                        protein_id = record.name
                        fasta = os.path.dirname(os.getcwd()) + f'{complete_genome_dir}/{phage[1]}/{phage[1]}.ffn'
                        hypothetical_fasta_str += f'>{protein_id} {product} \n{self.get_seq_in_fasta(protein_id, fasta)}\n'

                    if len(rbp_fasta_str) > 0:
                        file_name = f'{phage[1]}-rbp.ffn'

                        with open(os.path.join(rbp_prokka_dir, file_name), 'w') as rbp_fasta_file:
                            rbp_fasta_file.write(rbp_fasta_str)
                            
                    if len(hypothetical_fasta_str) > 0:
                        file_name = f'{phage[1]}-hypothetical.ffn'

                        with open(os.path.join(hypothetical_prokka_dir, file_name), 'w') as hypothetical_fasta_file:
                            hypothetical_fasta_file.write(hypothetical_fasta_str)
                            
            # Display progress
            ctr += 1
            self.__display_progress(ctr)