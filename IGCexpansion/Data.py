# A separate class to represent data structure
# Xiang Ji
# xji3@ncsu.edu

from Bio import SeqIO
import os, sys

class Data:
    def __init__(self, alignment_file, gene_to_orlg_file, two_sites = False):
        self.nsites             = None               # Number of sites in the alignment
        self.alignment_file     = alignment_file     # Multiple sequence alignment file location
        self.gene_to_orlg_file  = gene_to_orlg_file  # Gene ortholog mapping info file location

        self.name_to_seq        = None               # dictionary used to store sequence
        # name should be species name + gene name
        self.gene_to_orlg       = dict()             # dictionary used to store gene ortholog group info
        self.two_sites_name_to_seq = dict()            # For new two sites model

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
        if two_sites:
            for space in range(1, self.nsites):
                self.two_sites_name_to_seq[space] = self.get_two_sites_states(space)

    def is_alignment(self): # test if all sequences are of same length
        return len(set([len(self.name_to_seq[name]) for name in self.name_to_seq])) == 1
 

    def get_two_sites_states(self, space, data_type = 'nt'):
        # space = 1 represents neighboring sites
        if data_type == 'nt':
            obs_to_state = {'ACGT'[nt]:nt for nt in range(4)}
        else:
            sys.exit('The data_type is not supported in Data class.')

        assert(0 < space < self.nsites)
        new_name_to_pair_state = dict()
        for name in self.name_to_seq:
            seq = self.name_to_seq[name]
            new_name_to_pair_state[name] = [(obs_to_state[seq[i]], obs_to_state[seq[i + space]]) for i in range(len(seq) - space)]

        return new_name_to_pair_state
            
        


            

if __name__ == '__main__':
    gene_to_orlg_file = '../test/ADH1GeneToOrlg.txt'
    alignment_file = '../test/ADH1Alignment_test.txt'
    
    test = Data(alignment_file, gene_to_orlg_file)
    self = test
    print test.gene_to_orlg, test.name_to_seq, test.nsites
##    gene_to_orlg_file = '../test/Trigeneconv_GeneToOrlg.txt'
##    alignment_file = '../test/Trigeneconv_ADH_intron_input.fasta'
##    
##    test = Data(alignment_file, gene_to_orlg_file, True)
##    self = test
##    print test.gene_to_orlg
    
    gene_to_orlg_file = '../test/YDR418W_YEL054C_GeneToOrlg.txt'
    alignment_file = '../test/YDR418W_YEL054C_MG94_geo_10.0_Sim_8.fasta'
    
    test = Data(alignment_file, gene_to_orlg_file, True)
    self = test
    #print test.gene_to_orlg, test.name_to_seq
    name_to_pair_state = test.get_two_sites_states(0)

