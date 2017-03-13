# A separate file for Simulator
# This simulation should simulate multigene family evolution with
# point mutation and IGC(interlocus gene conversion) processes
# along a given species tree with gene duplication loss history and arbitrary family size
# This simulator only simulates cdna for now
# Xiang Ji
# xji3@ncsu.edu

from IGCTractModel import IGCTractModel
from Tree import Tree
from PMModel import PMModel
from Common import divide_configuration, draw_from_distribution
import numpy as np
import cPickle, os


class Simulator:

    def __init__(self, pm_model_name, x_pm, rate_variation,   # PM Model
                 x_IGC, pm_IGC,                               # IGC Tract Model
                 tree_newick, DupLosList,                     # Tree input
                 terminal_node_list, node_to_pos,             # Configuration input
                 gene_to_orlg_file, seq_file, log_file,       # output info
                 seed_file,                                   # random seed file
                 seq_index_file, nsites,                      # sequence index file, number of nt to simulate
                 cdna = True                                  # simulate cdna sequence only or not
                 ):

        self.tree = Tree(tree_newick, DupLosList, terminal_node_list, node_to_pos)
        self.node_to_pos = node_to_pos
        self.terminal_node_list = terminal_node_list                 
        self.root_by_dup = self.tree.is_duplication_node(self.tree.phylo_tree.root.name)
        self.PMModel = PMModel(pm_model_name, x_pm, rate_variation)
        self.IGCModel = IGCTractModel(x_IGC, [0, 1], pm_IGC)
        self.seq_index_file = seq_index_file # seq_index file location
        self.seq_index = None                # store sequence index information
        self.nsites    = nsites              # number of nucleotide to simulate, must agree with seq_index length
        self.cdna      = cdna

        self.seq_file  = seq_file            # output sequence file
        self.log_file  = log_file            # output log file
        self.seed_file = seed_file           # random seed file


        self.initiate()

    def __str__(self):  # overide for print function
        print self.PMModel
        print self.IGCModel
        print self.tree

        return 'IGC simulator output seq: ' + self.seq_file + '\n' + \
               'log file: ' + self.log_file + '\n'

    def initiate(self):
        self.get_seed_file()
        self.read_seq_index_file()
        

    def get_seed_file(self):
        if os.path.isfile(self.seed_file):
            prng = cPickle.load(open(self.seed_file, 'r'))
            np.random.set_state(prng.get_state())
        else:
            prng = np.random.RandomState()
            cPickle.dump(prng, open(self.seed_file, 'w+'))            

    def read_seq_index_file(self):
        # The index should have columns:
        # nt_index, codon #, codon site for coding sequence
        # nt_index,  -1,  -1 for noncoding sequence
        if self.seq_index_file == None:
            if self.cdna:
                assert(self.nsites %3 == 0)
                return np.array([[i + 1, int(floor((i + 0.5) / 3.0)) + 1, i + 1 - 3 * int(floor((i + 0.5) / 3.0))] for i in range(self.nsites)])
            else:
                return np.array([[i, -1, -1] for i in range(self.nsites)])
        else:
            seq_index = np.loadtxt(self.seq_index_file, dtype = int)
            if not self.cdna:
                seq_index[:, 1:] = -1
        self.seq_index = seq_index
        assert(len(self.seq_index) == self.nsites)

    def sim_root(self):
        if self.PMModel.name == 'HKY':
            distn = [self.PMModel.parameters['Pi_' + nt] for nt in 'ACGT']
            seq = draw_from_distribution(distn, self.nsites, 'ACGT')
        return seq
            
        

if __name__ == '__main__':
    gene_to_orlg_file = '../test/YDR418W_YEL054C_GeneToOrlg.txt'
    seq_file = '../test/YDR418W_YEL054C_Simulation.fasta'
    log_file = '../test/YDR418W_YEL054C_Simulation.log'
    seed_file = '../test/YDR418W_YEL054C_Simulation_seed.log'

    tree_newick = '../test/YeastTree.newick'
    DupLosList = '../test/YeastTestDupLost.txt'
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    seq_index_file = '../test/YDR418W_YEL054C_seq_index.txt'
    nsites = 489

    pm_model_name = 'HKY'
    x_pm = np.log([0.4, 0.5, 0.2, 9.2])
    rate_variation = False

    x_IGC = [2.0, 0.3]
    init_pm = 'One rate'
    tract_pm = 'One rate'
    pm_IGC = [init_pm, tract_pm]
    
    test = Simulator(pm_model_name, x_pm, rate_variation,
                     x_IGC, pm_IGC, tree_newick, DupLosList,
                     terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, log_file, seed_file, seq_index_file, nsites)

    self = test
    print test
    print test.sim_root()
