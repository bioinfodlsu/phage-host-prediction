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
        - min_len_keyword: Minimum length for a 
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
        Sets the 
        """
        self.inphared_gb = inphared_gb
        
    def set_inphared(self, inphared):
        self.inphared = inphared
           
    def set_no_cds_annot(self, no_cds_annot):
        self.no_cds_annot = no_cds_annot
        
    def set_annot_products(self, annot_products):
        self.annot_products = annot_products
        
    def set_rbp_products(self, rbp_products):
        self.rbp_products = rbp_products
        
    def set_hypothetical_proteins(self, hypothetical_proteins):
        self.hypothetical_proteins = hypothetical_proteins
        
    def set_token_delimiter(self, token_delimiter):
        self.token_delimiter = token_delimiter

    # ==================
    # Data Preprocessing
    # ==================
        
    def get_phages_unspec_host(self, inphared_unspec_host):
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
        unfiltered_hosts = set()
        for key, value in self.phages_unspec_host.items():
            for host in value[0]:
                unfiltered_hosts.add(host)
            for host in value[1]:
                unfiltered_hosts.add(host)
                
        self.unfiltered_hosts = unfiltered_hosts

        return unfiltered_hosts
    
    def get_genus_typo(self, filename):
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
        unfiltered_suspected_genera = set()

        for unfiltered_host in self.unfiltered_hosts:
            unfiltered_suspected_genera.add(self.__get_suspected_genus(candidate_regex, unfiltered_host))
            
        self.unfiltered_suspected_genera = unfiltered_suspected_genera

        return unfiltered_suspected_genera
    
    def get_excluded_hosts(self, *exclusion_files):
        excluded_hosts = set()
        
        for exclusion_file in exclusion_files:
            with open(exclusion_file, 'r') as file:
                for host in file:
                    excluded_hosts.add(host.rstrip('\n'))
        
        self.excluded_hosts = excluded_hosts
        
        return excluded_hosts
    
    def get_valid_hosts(self):
        self.valid_hosts = self.unfiltered_suspected_genera - self.excluded_hosts
        return self.valid_hosts
    
    def is_possible_misspelling(self, typo, keyword):
        return (len(keyword) >= self.min_len_keyword and 
                nltk.edit_distance(typo, keyword, transpositions = True) <= self.misspelling_threshold)
    
    def update_host_column(self, candidate_regex, inphared_unspec_host, inphared_augmented):
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
        unknown_aa_expr = r'[^ACDEFGHIKLMNPQRSTVWY]'
        if re.search(unknown_aa_expr, sequence):
            return True

        return False
 
    def __is_acceptable_length(self, sequence, lower_bound, upper_bound):
        return lower_bound <= len(sequence) and len(sequence) <= upper_bound
    
    def __is_rbp(self, product, sequence, lower_bound = -1, upper_bound = -1):
        if lower_bound == -1 and upper_bound == -1:
            return (product in self.rbp_products
                    and not self.__has_unknown_amino_acid(sequence))
        
        return (product in self.rbp_products
                and not self.__has_unknown_amino_acid(sequence)
                and self.__is_acceptable_length(sequence, lower_bound, upper_bound))
    
    def generate_rbp_len_distribution(self):
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
        return (product in self.hypothetical_proteins 
                and not self.__has_unknown_amino_acid(sequence) 
                and self.__is_acceptable_length(sequence, lower_bound, upper_bound))
    
    def generate_rbp_hypothetical_fasta(self, rbp_genbank_dir, hypothetical_genbank_dir, lower_bound, upper_bound):
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
        fasta = set()
        for file in os.listdir(fasta_dir):
            fasta.add(file[:-(1 + len(suffix) + len('.fasta'))])

        embeddings = set()
        for file in os.listdir(embed_dir):
            embeddings.add(file[:-(1 + len(suffix) + len('-embeddings.csv'))])

        return embeddings - fasta, fasta - embeddings
    
    def check_fasta_embeddings_per_protein(self, suffix, fasta_dir, embed_dir):
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
        for record in SeqIO.parse(fasta, 'fasta'):
            record_id = record.description.split(' ', 1)[0]
            
            if record_id == protein_id:
                return str(record.seq)
            
        return None
    
    def generate_rbp_hypothetical_nucleotide_prokka(self, complete_genome_dir, rbp_prokka_dir, hypothetical_prokka_dir,
                                                    lower_bound, upper_bound):
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