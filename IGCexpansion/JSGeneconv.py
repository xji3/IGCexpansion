# A separate file for JSLikelihood
# Xiang Ji
# xji3@ncsu.edu

from Data import Data
from Tree import Tree
from JSModel import JSModel
from Func import *
from copy import deepcopy

class JSGeneconv:
    def __init__(self, alignment_file, gene_to_orlg_file, # Data input
                 tree_newick, DupLosList,                 # Tree input
                 n_js, x_js, pm_model, n_orlg, IGC_pm,    # JSModel input
                 node_to_pos, terminal_node_list,         # Configuration input
                 root_by_dup = False):   # JSModel input
        self.tree = Tree(tree_newick, DupLosList)
        self.data = Data(alignment_file, gene_to_orlg_file)
        self.jsmodel = JSModel(n_js, x_js, pm_model, n_orlg, IGC_pm)
        self.node_to_pos = node_to_pos
        self.terminal_node_list = terminal_node_list
        self.root_by_dup = root_by_dup
        if self.root_by_dup:
            self.x = np.concatenate((np.log([0.1] * len(self.tree.edge_list)), x_js))
        else:
            self.x = np.concatenate((np.log([0.1] * (len(self.tree.edge_list) - 1)), x_js))

        self.tree.get_configurations(self.terminal_node_list, self.node_to_pos)
        assert(self.jsmodel.n_orlg == self.tree.n_orlg)  # make sure n_orlg got updated
        
    def unpack_x(self, x):
        self.x = x
        if self.root_by_dup:
            x_rate = x[:len(self.tree.edge_list)]
            x_js   = x[len(self.tree.edge_list):]
        else:
            x_rate = x[:(len(self.tree.edge_list) - 1)]
            x_js   = x[(len(self.tree.edge_list) - 1):]
        self.unpack_x_rates(x_rate)
        self.jsmodel.unpack_x_js(x_js)

    def unpack_x_rates(self, x_rate):
        if self.root_by_dup:
            self.tree.unpack_x_rates(x_rate)
        else: # add additional constrain on the brach from root to first duplication node
            # The first element in x_rate is the combined rate for both N0 related edges
            # First duplication node should be on one of the two N0 related edges
            # Constraint: N0_Outgroup = 9 * N0_D1
            total_rate = np.exp(x_rate[0])
            outgroup_rate = x_rate[0] + np.log(0.9)
            N0_D1_rate = x_rate[0] + np.log(0.1)
            translated_rate = []
            rate_iter = 1
            for edge in self.tree.edge_list:
                if 'N0' in edge:
                    if 'D1' in edge:
                        translated_rate.append(N0_D1_rate)
                    else:
                        translated_rate.append(outgroup_rate)
                else:
                    translated_rate.append(x_rate[rate_iter])
                    rate_iter += 1
            self.tree.unpack_x_rates(translated_rate)

    def get_scene(self):
        state_space_shape = self.jsmodel.state_space_shape
        process_definitions, conf_list = get_process_definitions(self.tree, self.jsmodel)
        self.tree.get_tree_process(conf_list)
        scene = dict(
            node_count = len(self.tree.node_to_num),
            process_count = len(process_definitions),
            state_space_shape = state_space_shape,
            tree = self.tree.tree_json,
            
            )

    
        



if __name__ == '__main__':
    gene_to_orlg_file = '../test/ADH1GeneToOrlg.txt'
    alignment_file = '../test/ADH1Alignment_test.txt'
    
    data = Data(alignment_file, gene_to_orlg_file)
    
    tree_newick = '../test/PrimateTest.newick'
    DupLosList = '../test/PrimateTestDupLost.txt'
    tree = Tree(tree_newick, DupLosList)

    node_to_pos = {'D1':0, 'D2':0, 'D3':1, 'D4':2, 'L1':2}
    terminal_node_list = ['Chinese_Tree_Shrew', 'Macaque', 'Olive_Baboon', 'Orangutan', 'Gorilla', 'Human']
    tree.get_configurations(terminal_node_list, node_to_pos)

    pm_model = 'HKY'
    x_js = np.log([0.3, 0.5, 0.2, 9.5, 4.9])
    n_orlg = 9
    IGC_pm = 'One rate'
    n_js = 5
    jsmodel = JSModel(n_js, x_js, pm_model, n_orlg, IGC_pm)

    test = JSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList, n_js, x_js, pm_model, n_orlg, IGC_pm,
                      node_to_pos, terminal_node_list)
    self = test
    x = np.concatenate((np.log([0.2] * (len(test.tree.edge_list) - 1)), x_js))
    test.unpack_x(x)
    process_definitions, conf_list = get_process_definitions(self.tree, self.jsmodel)
    
    conf_list = count_process(test.tree.node_to_conf)
    configuration = conf_list[0]
#    process = jsmodel.get_process_definition(configuration)
    
