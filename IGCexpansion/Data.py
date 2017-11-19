# A separate class to represent data structure
# Xiang Ji
# xji3@ncsu.edu

from Bio import SeqIO
import os, sys
from math import floor
from itertools import product
import numpy as np

class Data:
    def __init__(self, alignment_file, gene_to_orlg_file, seq_index_file = None,
                 two_sites = False, space_list = None, cdna = False, allow_same_codon = False):
        self.nsites             = None               # Number of sites in the alignment
        self.alignment_file     = alignment_file     # Multiple sequence alignment file location
        self.gene_to_orlg_file  = gene_to_orlg_file  # Gene ortholog mapping info file location
        self.space_list         = space_list

        self.name_to_seq        = None               # dictionary used to store sequence
        # name should be species name + gene name
        self.gene_to_orlg       = dict()             # dictionary used to store gene ortholog group info
        self.two_sites_name_to_seq = dict()          # For new two sites model
        self.max_space          = None               # space calculation control, its value is defined in read_seq_index_file function
        self.seq_index          = None
        self.seq_index_file     = seq_index_file
        # An index of the positions in the MSA (for cases where indels are removed or unconsidered)
        self.idx_to_pos         = dict()             # {self.seq_index[i, 0]:i}
        self.space_idx_pairs    = dict()
        self.cdna               = cdna               # a bool indicator for coding/non-coding sequences
        self.codon_site_to_idx  = None               # a dictionary to store idx of subsequence of different codon site
        self.codon_site_to_pos  = None               # a dictionary to store pos of subsequence of different codon site
        self.allow_same_codon   = allow_same_codon   # a bool indicator for coding sequence when used in pair site models
        self.space_codon_site_pair = None            # a dictionary to store available codon pair for given space. This is my lazy compromise.
        

        self.get_gene_to_orlg()
        self.get_data(two_sites)
        assert(self.self_check())

    def get_gene_to_orlg(self):
        assert(os.path.isfile(self.gene_to_orlg_file))
        with open(self.gene_to_orlg_file, 'rb') as f:
            for line in f:
                items = line.split()
                if items:
                    gene = items[0]
                    orlg = int(items[1])
                    self.gene_to_orlg[gene] = orlg

    def self_check(self):
        check_status = True
        if self.allow_same_codon:
            check_status = check_status and self.cdna
        max_space = max(self.seq_index[:, 0]) - min(self.seq_index[:, 0]) + 1
        check_status = check_status and all([0 < space < self.max_space + 1 for space in self.space_list])

        return check_status

        
    def read_seq_index_file(self, seq_index_file):
        # The index should have columns:
        # nt_index, codon #, codon site for coding sequence
        # nt_index,  -1,  -1 for noncoding sequence
        if seq_index_file == None:
            self.max_space = self.nsites - 1
            if self.cdna:
                assert(self.nsites %3 == 0)
                return np.array([[i + 1, int(floor((i + 0.5) / 3.0)) + 1, i + 1 - 3 * int(floor((i + 0.5) / 3.0))] for i in range(self.nsites)])
            else:
                return np.array([[i, -1, -1] for i in range(self.nsites)])
        else:
            seq_index = np.loadtxt(seq_index_file, dtype = int)
            self.max_space = max(seq_index[:, 0]) - min(seq_index[:, 0])
            if not self.cdna:
                seq_index[:, 1:] = -1
            return seq_index

        
    def get_data(self, two_sites):
        assert(os.path.isfile(self.alignment_file))
        seq_dict = SeqIO.to_dict(SeqIO.parse( self.alignment_file, "fasta" ))
        self.name_to_seq = {name:str(seq_dict[name].seq) for name in seq_dict.keys()}

        assert(self.is_alignment)
        self.nsites = len(self.name_to_seq[self.name_to_seq.keys()[0]])

        # Now assign self.seq_index
        self.seq_index = self.read_seq_index_file(self.seq_index_file)
        
        assert(len(self.seq_index) == self.nsites)
        self.idx_to_pos = {self.seq_index[i, 0]:i for i in range(self.seq_index.shape[0])}
        # and the seq_index has to be increasing order
        assert(all(earlier <= later for earlier, later in zip(self.seq_index[:, 0], self.seq_index[1:, 0])))

        if self.cdna:
            self.sort_codon_site()

        if self.space_list == None:
            self.space_list = self.get_possible_space_list()
        
        if two_sites:
            self.get_two_sites_states()

            

    def is_alignment(self): # test if all sequences are of same length
        return len(set([len(self.name_to_seq[name]) for name in self.name_to_seq])) == 1

    def sort_codon_site(self): # for cdna sequences divide the sequence into three subsequences
        assert(self.cdna) # this function should only be used for cdna sequences
        self.codon_site_to_idx = {1:[], 2:[], 3:[]}
        self.codon_site_to_pos = {1:[], 2:[], 3:[]}
        for i in range(self.seq_index.shape[0]):
            nt_idx, codon_idx, codon_site = self.seq_index[i, :]
            self.codon_site_to_idx[codon_site].append(nt_idx)
            self.codon_site_to_pos[codon_site].append(i)
            
        

    def get_possible_space_list(self):
        possible_space_list = list()
        if self.cdna:
            for i, j in product(range(1, 4), repeat = 2):
                self.space_idx_pairs[(i, j)] = dict()
                
            if self.allow_same_codon:
                for i in range(self.seq_index.shape[0] - 1):
                    for j in range(i + 1, self.seq_index.shape[0]):
                        nt_idx_i, codon_idx_i, codon_site_i = self.seq_index[i, :]
                        nt_idx_j, codon_idx_j, codon_site_j = self.seq_index[j, :]
                        space = nt_idx_j - nt_idx_i
                        if space in self.space_idx_pairs[(codon_site_i, codon_site_j)]:
                            self.space_idx_pairs[(codon_site_i, codon_site_j)][space].append((i, j))
                        else:
                            self.space_idx_pairs[(codon_site_i, codon_site_j)][space] = [(i, j)]
                        possible_space_list.append(space)
            else:
                for i in range(self.seq_index.shape[0] - 1):
                    for j in range(i + 1, self.seq_index.shape[0]):
                        nt_idx_i, codon_idx_i, codon_site_i = self.seq_index[i, :]
                        nt_idx_j, codon_idx_j, codon_site_j = self.seq_index[j, :]
                        space = nt_idx_j - nt_idx_i
                        if codon_idx_i != codon_idx_j:
                            if space in self.space_idx_pairs[(codon_site_i, codon_site_j)]:
                                self.space_idx_pairs[(codon_site_i, codon_site_j)][space].append((i, j))
                            else:
                                self.space_idx_pairs[(codon_site_i, codon_site_j)][space] = [(i, j)]
                            possible_space_list.append(self.seq_index[j,0] - self.seq_index[i, 0])
        else:
            for i in range(self.seq_index.shape[0] - 1):
                for j in range(i + 1, self.seq_index.shape[0]):
                    space = self.seq_index[j,0] - self.seq_index[i, 0]
                    if space in self.space_idx_pairs:
                        self.space_idx_pairs[space].append((i, j))
                    else:
                        self.space_idx_pairs[space] = [(i, j)]
                    possible_space_list.append(space)
            

        possible_space_list = list(set(possible_space_list))

        # My lazy compromise
        if self.cdna:
            self.space_codon_site_pair = {n:[codon_site_pair for codon_site_pair in product(range(1, 4), repeat = 2) if n in self.space_idx_pairs[codon_site_pair]] for n in possible_space_list}
        return possible_space_list

            
 

    def get_two_sites_states(self, data_type = 'nt'):
        # space = 1 represents neighboring sites
        if data_type == 'nt':
            obs_to_state = {'ACGT'[nt]:nt for nt in range(4)}
            obs_to_state['-'] = -1
            obs_to_state['N'] = -1
        else:
            sys.exit('The data_type is not supported in Data class.')

        if self.cdna:
            for i, j in product(range(1, 4), repeat = 2):
                self.two_sites_name_to_seq[(i, j)] = dict()
                for space in self.space_idx_pairs[(i, j)]:
                    new_name_to_pair_state = dict()
                    for name in self.name_to_seq:
                        seq = self.name_to_seq[name]
                        ps_state_list = [(obs_to_state[seq[idx_pair[0]]], obs_to_state[seq[idx_pair[1]]]) for idx_pair in self.space_idx_pairs[(i, j)][space]]
                        new_name_to_pair_state[name] = ps_state_list            
                    self.two_sites_name_to_seq[(i, j)][space] = new_name_to_pair_state
        else:
            for space in self.space_list:
                new_name_to_pair_state = dict()
                for name in self.name_to_seq:
                    seq = self.name_to_seq[name]
                    ps_state_list = [(obs_to_state[seq[idx_pair[0]]], obs_to_state[seq[idx_pair[1]]]) for idx_pair in self.space_idx_pairs[space]]
                    new_name_to_pair_state[name] = ps_state_list
                self.two_sites_name_to_seq[space] = new_name_to_pair_state

        


            

if __name__ == '__main__':
##    gene_to_orlg_file = '../test/ADH1GeneToOrlg.txt'
##    alignment_file = '../test/ADH1Alignment_test.txt'
##    
##    test = Data(alignment_file, gene_to_orlg_file)
##    self = test
##    print test.gene_to_orlg, test.name_to_seq, test.nsites
##    gene_to_orlg_file = '../test/Trigeneconv_GeneToOrlg.txt'
##    alignment_file = '../test/Trigeneconv_ADH_intron_input.fasta'
##    
##    test = Data(alignment_file, gene_to_orlg_file, True)
##    self = test
##    print test.gene_to_orlg
    
    gene_to_orlg_file = '../test/YDR418W_YEL054C_GeneToOrlg.txt'
    alignment_file = '../test/YDR418W_YEL054C_MG94_geo_10.0_Sim_8.fasta'
    seq_index_file = '../test/YDR418W_YEL054C_seq_index.txt'
    seq_index_file = None
    space_list = None
    
    test = Data(alignment_file, gene_to_orlg_file, two_sites = True,
                space_list = space_list, seq_index_file = seq_index_file, cdna = True, allow_same_codon = True)
    print test.seq_index
    self = test
    #print test.gene_to_orlg, test.name_to_seq
    #name_to_pair_state = test.get_two_sites_states(1)
    possible_space_list = test.get_possible_space_list()
    print len(possible_space_list), test.nsites

