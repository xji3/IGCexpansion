# A separate file for JSLikelihood
# Xiang Ji
# xji3@ncsu.edu

from Data import Data
from Tree import Tree
from JSModel import JSModel
from Func import *

class JSGeneconv:
    def __init__(self, alignment_file, gene_to_orlg_file, # Data input
                 tree_newick, DupLosList,                 # Tree input
                 n_js, x_js, pm_model, n_orlg, IGC_pm):   # JSModel input
        self.tree = Tree(tree_newick, DupLosList)
        self.data = Data(alignment_file, gnee_to_orlg_file)
        self.jsmodel = JSModel(n_js, x_js, pm_model, n_orlg, IGC_pm)



if __name__ == '__main__':
    print
