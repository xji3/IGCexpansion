# A separate file for general functions
# Xiang Ji
# xji3@ncsu.edu

from Data import Data
from Tree import Tree
from JSModel import JSModel
import numpy as np
from copy import deepcopy


def get_iid_observations(data, tree, nsites, data_type = 'nt'):
    assert(isinstance(data, Data) and isinstance(tree, Tree))
    # Generate iid_observations from data class and tree class
    if data_type == 'nt':
        obs_to_state = {nt:'ACGT'.index(nt) for nt in 'ACGT'}

    observable_names = data.gene_to_orlg.keys()  # don't have to use all sequences
    name_node_axes = []
    visited_terminal_nodes = []
    for name in observable_names:
        terminal_node, gene_name = name.split('__')
        orlg = data.gene_to_orlg[name]
        orlg_to_pos = tree.divide_configuration(tree.node_to_conf[terminal_node])
        assert(orlg in orlg_to_pos['extent'])
        for pos in orlg_to_pos['extent'][orlg]:
            name_node_axes.append([name, tree.node_to_num[terminal_node], pos])

        if not terminal_node in visited_terminal_nodes:
            for pos in orlg_to_pos['distinct']:
                name_node_axes.append(['distinct', tree.node_to_num[terminal_node], pos])
            visited_terminal_nodes.append(terminal_node)
        
    name_node_axes = sorted(name_node_axes, key = lambda name:name[1]) # sort by terminal_node
    observable_nodes = [item[1] for item in name_node_axes]
    observable_axes  = [item[2] for item in name_node_axes]
    observed_names   = [item[0] for item in name_node_axes]

    iid_observations = []
    for site in range(nsites):
        observations = []
        for name in observed_names:
            if name == 'distinct':
                observation = -1
            elif data.name_to_seq[name][site] == 'n' or data.name_to_seq[name][site] == '-':
                observation = -1
            else:
                observation = obs_to_state[data.name_to_seq[name][site].upper()]
            observations.append(observation)
        iid_observations.append(observations)

    return observable_nodes, observable_axes, iid_observations


def count_process(node_to_conf):
    conf_list = []
    for node in node_to_conf:
        if node_to_conf[node] in conf_list:
            continue
        else:
            conf_list.append(deepcopy(node_to_conf[node]))

    return conf_list
    
def get_process_definitions(tree, jsmodel, proportions = False):
    assert(isinstance(tree, Tree) and isinstance(jsmodel, JSModel))
    conf_list = count_process(tree.node_to_conf)
    process_definitions = []
    for conf in conf_list:
        #print conf
        process = jsmodel.get_process_definition(conf, proportions)
        process_definitions.append(process)
    return process_definitions, conf_list

def get_directional_process_definitions(tree, jsmodel, orlg_pair):
    assert(isinstance(tree, Tree) and isinstance(jsmodel, JSModel))
    conf_list = count_process(tree.node_to_conf)
    process_definitions = []
    for conf in conf_list:
        #print conf
        process = jsmodel.get_directional_process_definition(conf, orlg_pair)
        process_definitions.append(process)
    return process_definitions, conf_list


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
    print tree.dup_events
    for i in tree.node_to_conf:
        if i in terminal_node_list:
            print i, tree.node_to_conf[i]

    pm_model = 'HKY'
    x_js = np.log([0.3, 0.5, 0.2, 9.5, 4.9])
    n_orlg = tree.n_orlg
    IGC_pm = 'One rate'
    n_js = 9
    jsmodel = JSModel(n_js, x_js, pm_model, n_orlg, IGC_pm)
    



# function get_iid_observations()
    data_type = 'nt'
    nsites = 5
    
    observable_nodes, observable_axes, iid_observations = get_iid_observations(data, tree, nsites, data_type)

        
# function get_process_definitions
    process_definitions, conf_list = get_process_definitions(tree, jsmodel)

# test tree.get_tree_process(conf_list):
    tree.get_tree_process(conf_list)
    print tree.tree_json
