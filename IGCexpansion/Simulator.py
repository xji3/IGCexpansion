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
import numpy as np


class Simulator:

    def __init__(self, pm_model_name, x_pm, rate_variation,   # PM Model
                 x_IGC, pm_IGC,                               # IGC Tract Model
                 tree_newick, DupLosList,                     # Tree input
                 terminal_node_list, node_to_pos,             # Configuration input
                 gene_to_orlg_file, seq_file, log_file,       # output info
                 seq_index_file                               # sequence index file
                 ):

        self.tree = Tree(tree_newick, DupLosList, terminal_node_list, node_to_pos)
        self.node_to_pos = node_to_pos
        self.terminal_node_list = terminal_node_list                 
        self.root_by_dup = self.tree.is_duplication_node(self.tree.phylo_tree.root.name)
        self.PMModel = PMModel(pm_model_name, x_pm, rate_variation)
        self.IGCModel = IGCTractModel(x_IGC, [0, 1], pm_IGC)

        self.seq_file  = seq_file            # output sequence file
        self.log_file  = log_file            # output log file

    def __str__(self):  # overide for print function
        print self.PMModel
        print self.IGCModel
        print self.tree

        return 'IGC simulator output seq: ' + self.seq_file + '\n' + \
               'log file: ' + self.log_file + '\n'



if __name__ == '__main__':
    gene_to_orlg_file = '../test/YDR418W_YEL054C_GeneToOrlg.txt'
    seq_file = '../test/YDR418W_YEL054C_Simulation.fasta'
    log_file = '../test/YDR418W_YEL054C_Simulation.log'

    tree_newick = '../test/YeastTree.newick'
    DupLosList = '../test/YeastTestDupLost.txt'
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    seq_index_file = '../test/YDR418W_YEL054C_seq_index.txt'

    pm_model_name = 'HKY'
    x_pm = np.log([0.4, 0.5, 0.2, 9.2])
    rate_variation = False

    x_IGC = [2.0, 0.3]
    init_pm = 'One rate'
    tract_pm = 'One rate'
    pm_IGC = [init_pm, tract_pm]
    
    test = Simulator(pm_model_name, x_pm, rate_variation,
                     x_IGC, pm_IGC, tree_newick, DupLosList,
                     terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, log_file, seq_index_file)

    self = test
    print test
