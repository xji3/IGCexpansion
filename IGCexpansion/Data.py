# A separate class to represent data structure
# Xiang Ji
# xji3@ncsu.edu

from Bio import SeqIO
import os, sys

class Data:
    def __init__(self, alignment_file, gene_to_orlg_file, seq_index = None, two_sites = False, space_list = None):
        self.nsites             = None               # Number of sites in the alignment
        self.alignment_file     = alignment_file     # Multiple sequence alignment file location
        self.gene_to_orlg_file  = gene_to_orlg_file  # Gene ortholog mapping info file location
        self.space_list         = space_list

        self.name_to_seq        = None               # dictionary used to store sequence
        # name should be species name + gene name
        self.gene_to_orlg       = dict()             # dictionary used to store gene ortholog group info
        self.two_sites_name_to_seq = dict()          # For new two sites model
        self.seq_index          = seq_index          # An index of the positions in the MSA (for cases where indels are removed or unconsidered)

        self.get_gene_to_orlg()
        self.get_data(two_sites)

    def get_gene_to_orlg(self):
        assert(os.path.isfile(self.gene_to_orlg_file))
        with open(self.gene_to_orlg_file, 'rb') as f:
            for line in f:
                items = line.split()
                if items:
                    gene = items[0]
                    orlg = int(items[1])
                    self.gene_to_orlg[gene] = orlg
 
    def get_data(self, two_sites):
        assert(os.path.isfile(self.alignment_file))
        seq_dict = SeqIO.to_dict(SeqIO.parse( self.alignment_file, "fasta" ))
        self.name_to_seq = {name:str(seq_dict[name].seq) for name in seq_dict.keys()}

        assert(self.is_alignment)
        self.nsites = len(self.name_to_seq[self.name_to_seq.keys()[0]])
        if self.seq_index == None:
            self.seq_index = range(self.nsites)
        assert(len(self.seq_index) == self.nsites)
        # and the seq_index has to be increasing order
        assert(all(earlier <= later for earlier, later in zip(self.seq_index, self.seq_index[1:])))
        
        if two_sites:
            if self.space_list == None:
                self.space_list = self.get_possible_space_list()
            for space in self.space_list:
                self.two_sites_name_to_seq[space] = self.get_two_sites_states(space)

    def is_alignment(self): # test if all sequences are of same length
        return len(set([len(self.name_to_seq[name]) for name in self.name_to_seq])) == 1

    def get_possible_space_list(self):
        possible_space_list = list()
        for i in range(len(self.seq_index) - 1):
            for j in range(i + 1, len(self.seq_index)):
                possible_space_list.append(self.seq_index[j] - self.seq_index[i])

        possible_space_list = list(set(possible_space_list))
        return possible_space_list
 

    def get_two_sites_states(self, space, data_type = 'nt'):
        # space = 1 represents neighboring sites
        if data_type == 'nt':
            obs_to_state = {'ACGT'[nt]:nt for nt in range(4)}
        else:
            sys.exit('The data_type is not supported in Data class.')

        if not 0 < space < self.nsites:
            sys.exit('Change space please. Minimum = 1, Maximum = ' + str(self.nsites))
        new_name_to_pair_state = dict()
        for name in self.name_to_seq:
            seq = self.name_to_seq[name]
            ps_state_list = []
            for idx in self.seq_index:
                if idx + space in self.seq_index:
                    ps_state_list.append((obs_to_state[seq[idx]], obs_to_state[seq[idx + space]]))
            new_name_to_pair_state[name] = ps_state_list

        return new_name_to_pair_state
            
        


            

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
    seq_index = range(489)
    space_list = None
    
    test = Data(alignment_file, gene_to_orlg_file, two_sites = True, space_list = space_list, seq_index = seq_index)
    self = test
    #print test.gene_to_orlg, test.name_to_seq
    name_to_pair_state = test.get_two_sites_states(1)
    possible_space_list = test.get_possible_space_list()
    print len(possible_space_list), test.nsites

