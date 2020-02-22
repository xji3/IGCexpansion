# A separate file for general functions
# Xiang Ji
# xji3@ncsu.edu

from Data import Data
from Tree import Tree
from JSModel import JSModel
import numpy as np
from copy import deepcopy
from itertools import product
from Common import *


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
        
    #name_node_axes = sorted(name_node_axes, key = lambda name:name[1]) # sort by terminal_node
    name_node_axes = sorted(name_node_axes) 
    observable_nodes = [item[1] for item in name_node_axes]
    observable_axes  = [item[2] for item in name_node_axes]
    observed_names   = [item[0] for item in name_node_axes]

    if data.cdna:
        iid_observations = {i:[] for i in range(1, 4)}
        for codon_site in range(1, 4):
            if nsites == None:
                iter_max = len(data.codon_site_to_pos[codon_site])
            else:
                iter_max = nsites
            for site in range(iter_max):
                observations = []
                for name in observed_names:
                    observed_nt = data.name_to_seq[name][data.codon_site_to_pos[codon_site][site]]
                    if name == 'distinct':
                        observation = -1
                    elif observed_nt.upper() == 'N' or observed_nt == '-':
                        observation = -1
                    else:
                        observation = obs_to_state[observed_nt.upper()]
                    observations.append(observation)
                
                iid_observations[codon_site].append(observations)
    else:
        assert(nsites)
        iid_observations = []
        for site in range(nsites):
            observations = []
            for name in observed_names:
                if name == 'distinct':
                    observation = -1
                elif data.name_to_seq[name][site].upper() == 'N' or data.name_to_seq[name][site] == '-':
                    observation = -1
                else:
                    observation = obs_to_state[data.name_to_seq[name][site].upper()]
                observations.append(observation)
            
            iid_observations.append(observations)

    return observable_nodes, observable_axes, iid_observations


def get_PS_iid_observations(data, tree, nsites, n, codon_site_pair = None, data_type = 'nt'):
    assert(isinstance(data, Data) and isinstance(tree, Tree))
    if codon_site_pair == None:
        assert(not data.cdna)
        assert(nsites  <= len(data.space_idx_pairs[n]))
    else:
        assert(type(codon_site_pair) == tuple and len(codon_site_pair) == 2 and all([ 0 < i < 4 for i in codon_site_pair]))
        assert(n in data.two_sites_name_to_seq[codon_site_pair])
        assert(nsites  <= len(data.space_idx_pairs[codon_site_pair][n]))
    
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
            name_node_axes.append([name, tree.node_to_num[terminal_node], pos * 2])  # pos*2 for (4, 4, 4, 4), pos for (16, 16)
            name_node_axes.append([name, tree.node_to_num[terminal_node], pos * 2 + 1])

        if not terminal_node in visited_terminal_nodes:
            for pos in orlg_to_pos['distinct']:
                name_node_axes.append(['distinct', tree.node_to_num[terminal_node], pos])
            visited_terminal_nodes.append(terminal_node)
        
    name_node_axes = sorted(name_node_axes, key = lambda name:name[1]) # sort by terminal_node
    observable_nodes = [item[1] for item in name_node_axes]
    observable_axes  = [item[2] for item in name_node_axes]
    observed_names   = [item[0] for item in name_node_axes]

    # lazy compromise for changing state space shape from (16,16) back to (4,4,4,4)
    lazy_observed_names = observed_names[0::2]

    iid_observations = []
    for site in range(nsites):
        observations = []
        for name in lazy_observed_names:
            if name == 'distinct':
                sys.exit('There should not be distinction!')
            else:
                if codon_site_pair == None:
                    observation = data.two_sites_name_to_seq[n][name][site]#translate_two_nt_to_one_state(data.two_sites_name_to_seq[n][name][site])
                else:
                    observation = data.two_sites_name_to_seq[codon_site_pair][n][name][site]#translate_two_nt_to_one_state(data.two_sites_name_to_seq[codon_site_pair][n][name][site])
            #observations.append(observation) # this is for (16, 16)
            observations.extend(observation)  # this is for (4,4,4,4)
        #print observations, observed_names
        iid_observations.append(observations)

    return observable_nodes, observable_axes, iid_observations

def get_all_PS_iid_observations(data, tree, data_type = 'nt', nsites = None):
    assert(data_type == 'nt')
    if data.cdna:
        iid_observations = {(i, j):dict() for i, j in product(range(1, 4), repeat = 2)}
        for codon_site_pair in iid_observations:
            for n in data.two_sites_name_to_seq[codon_site_pair]:
                if nsites == None:
                    used_nsites = len(data.space_idx_pairs[codon_site_pair][n])
                else:
                    used_nsites = min([nsites, len(data.space_idx_pairs[codon_site_pair][n])])
                iid_observable_nodes, observable_axes, single_iid_observations = get_PS_iid_observations(data, tree, used_nsites, n, codon_site_pair = codon_site_pair, data_type = data_type)
                iid_observations[codon_site_pair][n] = single_iid_observations
    else:
        iid_observations = {n:[] for n in data.space_list}
        for n in data.space_list:
            if nsites == None:
                used_nsites = len(data.space_idx_pairs[n])
            else:
                used_nsites = min([nsites, len(data.space_idx_pairs[n])])
            iid_observable_nodes, observable_axes, single_iid_observations = get_PS_iid_observations(data, tree, used_nsites, n, codon_site_pair = None, data_type = data_type)
            iid_observations[n] = single_iid_observations
    return iid_observable_nodes, observable_axes, iid_observations

def count_process(node_to_conf):
    conf_list = []
    for node in node_to_conf:
        if node_to_conf[node] in conf_list:
            continue
        else:
            conf_list.append(deepcopy(node_to_conf[node]))

    return sorted(conf_list)
    
def get_process_definitions(tree, jsmodel, proportions = False, codon_site = 1):
    assert(isinstance(tree, Tree) and isinstance(jsmodel, JSModel))
    conf_list = count_process(tree.node_to_conf)
    process_definitions = []
    for conf in conf_list:
        #print conf
        process = jsmodel.get_process_definition(conf, proportions, codon_site)
        process_definitions.append(process)
    return process_definitions, conf_list

def get_mutation_reduction_definitions(tree, jsmodel, codon_site = 1):
    assert(isinstance(tree, Tree) and isinstance(jsmodel, JSModel))
    conf_list = count_process(tree.node_to_conf)
    process_definitions = []
    for conf in conf_list:
        #print conf
        process = jsmodel.get_mutation_reduction_definition(conf, codon_site)
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
    alignment_file = '../test/Concatenated_all_exon.fasta'
    
    
    
    tree_newick = '../test/PrimateTest.newick'
    DupLosList = '../test/PrimateTestDupLost.txt'


    node_to_pos = {'D1':0, 'D2':0, 'D3':1, 'D4':2, 'L1':2}
    terminal_node_list = ['Chinese_Tree_Shrew', 'Macaque', 'Olive_Baboon', 'Orangutan', 'Gorilla', 'Human']
    tree = Tree(tree_newick, DupLosList, terminal_node_list, node_to_pos)
    data = Data(alignment_file, gene_to_orlg_file)
    
    print (tree.dup_events)
    for i in tree.node_to_conf:
        if i in terminal_node_list:
            print(i, tree.node_to_conf[i])

    pm_model = 'HKY'
    x_js = np.log([0.3, 0.5, 0.2, 9.5, 4.9])
    n_orlg = tree.n_orlg
    IGC_pm = 'One rate'
    n_js = tree.n_js
    jsmodel = JSModel(n_js, x_js, pm_model, n_orlg, IGC_pm)

    conf_list = count_process(tree.node_to_conf)
    accessible_orlg_pair = get_accessible_orlg_pair(conf_list)
    print(accessible_orlg_pair, len(accessible_orlg_pair))
    



### function get_iid_observations()
##    data_type = 'nt'
##    nsites = 5
##    
##    observable_nodes, observable_axes, iid_observations = get_iid_observations(data, tree, nsites, data_type)
##
##        
### function get_process_definitions
##    process_definitions, conf_list = get_process_definitions(tree, jsmodel)
##
### test tree.get_tree_process(conf_list):
##    tree.get_tree_process(conf_list)
##    print tree.tree_json
