# A separate class to represent data structure
# Xiang Ji
# xji3@ncsu.edu

from Bio import SeqIO
import os

class Data:
    def __init__(self, alignment_file, gene_to_orlg_file):
        self.nsites             = None               # Number of sites in the alignment
        self.alignment_file     = alignment_file     # Multiple sequence alignment file location
        self.gene_to_orlg_file  = gene_to_orlg_file  # Gene ortholog mapping info file location

        self.name_to_seq        = None               # dictionary used to store sequence
        # name should be species name + gene name
        self.gene_to_orlg       = dict()             # dictionary used to store gene ortholog group info

        self.get_gene_to_orlg()
        self.get_data()

    def get_gene_to_orlg(self):
        assert(os.path.isfile(self.gene_to_orlg_file))
        with open(self.gene_to_orlg_file, 'rb') as f:
            for line in f:
                items = line.split()
                if items:
                    gene = items[0]
                    orlg = int(items[1])
                    self.gene_to_orlg[gene] = orlg
 
    def get_data(self):
        assert(os.path.isfile(self.alignment_file))
        seq_dict = SeqIO.to_dict(SeqIO.parse( self.alignment_file, "fasta" ))
        self.name_to_seq = {name:str(seq_dict[name].seq) for name in seq_dict.keys()}

        assert(self.is_alignment)
        self.nsites = len(self.name_to_seq[self.name_to_seq.keys()[0]])

    def is_alignment(self): # test if all sequences are of same length
        return len(set([len(self.name_to_seq[name]) for name in self.name_to_seq])) == 1
 


if __name__ == '__main__':
    gene_to_orlg_file = '../test/ADH1GeneToOrlg.txt'
    alignment_file = '../test/ADH1Alignment_test.txt'
    
    test = Data(alignment_file, gene_to_orlg_file)
    self = test
    print test.gene_to_orlg, test.name_to_seq
