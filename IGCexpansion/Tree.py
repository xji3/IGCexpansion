# A separate class to represent tree structure
# Xiang Ji
# xji3@ncsu.edu
from operator import itemgetter
from itertools import groupby
from Bio import Phylo
import networkx as nx
import os, sys
import numpy as np
from copy import deepcopy
from Common import *

class Tree:
    def __init__(self, tree_newick, DupLosList, terminal_node_list, node_to_pos):
        self.newicktree         = tree_newick   # newick tree file location
        self.duploslist         = DupLosList    # duplication loss nodes file location
        # Tree topology related variable
        # The tree for this project is fixed, but keep the read-in feature anyway
        self.phylo_tree         = None        # Phylo tree structure
        self.tree_json          = None        # store the tree dictionary used for json likelihood package parsing
        self.edge_to_blen       = None        # dictionary store the unpacked tree branch length information {(node_from, node_to):blen}
        self.edge_list          = None        # kept all edges in the same order with x_rates
        self.node_to_num        = None        # dictionary used for translating tree info from self.edge_to_blen to self.tree
        self.num_to_node        = None        # dictionary used for translating tree info from self.tree to self.edge_to_blen
        self.node_to_dup        = dict()      # used to keep which orlg numbers are new duplicates on each duplication node

        # info for configurations
        self.terminal_node_list = terminal_node_list
        self.node_to_pos        = node_to_pos
        self.node_to_conf       = dict()      # A dictionary store configurations on each node

        # Speciation node starts with N, Duplication node with D, Loss node with L
        self.dup_events         = dict()      # A dictionary stores duplication events: ortholog group in key gives birth to the two ortholog groups in the content list
        
        self.n_orlg             = 0           # number of ortholog groups
        self.visited_DL_nodes   = list()
        self.n_js               = None
        self.get_tree()



    def get_tree(self):
        tree = Phylo.read( self.newicktree, "newick")
        self.phylo_tree = tree.as_phyloxml(rooted = 'True')
        self.add_duplos_nodes()
        #set node number for nonterminal nodes and specify root node

        self.get_tree_json()
        # get process function is implemented in Func.py

        # get node_to_conf
        self.get_configurations()

        # if 'Root' node is added as root node
        # there must be a duplication node directly after it and the Root node has no outgroup
        # delete the Root node after all configuration
        if self.phylo_tree.root.name == 'Root':
            assert(len(self.phylo_tree.root.clades) == 1)
            print 'Now remove root node and start with first duplication node'
            self.node_to_conf.pop(self.phylo_tree.root.name)
            self.root_with_duplication(self.phylo_tree.root.clades[0].name)
            
        # Warning:
        # cannot remove root node if it has multiple paralogs


    def get_tree_json(self):
        tree_nx = Phylo.to_networkx(self.phylo_tree)

        triples = [(u.name, v.name, d['weight']) for (u, v, d) in tree_nx.edges(data = True)] # data = True to have the blen as 'weight'

        T = nx.DiGraph()
        edge_to_blen = {}
        for va, vb, blen in triples:
            edge = (va, vb)
            T.add_edge(*edge)
            edge_to_blen[edge] = blen

        self.edge_to_blen = edge_to_blen

        # Now sort all nodes according to the degree where degree is as defined in graph theory
        all_nodes = sorted(T.degree().items(), key = lambda node: node[1])
        self.node_to_num = {n[0]:i for i, n in enumerate(all_nodes)}
        self.num_to_node = {i:n[0] for i, n in enumerate(all_nodes)}
        edge_list = sorted(edge_to_blen.keys())
         
        # Now setup self.tree dictionary
        tree_row = [self.node_to_num[na] for na, nb in edge_list]
        tree_col = [self.node_to_num[nb] for na, nb in edge_list]
        self.edge_list = edge_list

        self.tree_json = dict(
            row_nodes = tree_row,
            column_nodes = tree_col,
            #process = tree_process,
            edge_rate_scaling_factors = np.ones(len(tree_row))
            )

        
    def get_tree_process(self, conf_list):
        tree_process = []
        for edge in self.edge_list:
            parent_node, child_node = edge
            conf = self.node_to_conf[parent_node]
            tree_process.append(conf_list.index(conf))
        self.tree_json['edge_processes'] = tree_process

            
    def add_duplos_nodes(self):
        assert(os.path.isfile(self.duploslist))
        with open(self.duploslist, 'rb') as f:
            for line in f:
                items = line.split()
                if items:
                    branch = items[0]
                    first_ = branch.find('_')
                    father_node_name = branch[:first_]
                    child_node_name = branch[(first_ + 1):]
                    if father_node_name == 'Root' and self.phylo_tree.root.name != 'Root':
                        self.add_root()
                    
                    for add_node in items[1:]:
                        self.add_node(father_node_name, child_node_name, add_node)
                        father_node_name = add_node


            
                    
    def add_node(self, father_node_name, child_node_name, add_node):
        child_clade = self.find_clade(child_node_name)

        father_clade = self.find_parent_clade(child_node_name)
        assert(father_clade.name == father_node_name)

        new_clade = Phylo.BaseTree.Clade(name = add_node, clades = [child_clade])
        father_clade.clades.remove(child_clade)
        father_clade.clades.append(new_clade)

    def add_root(self, add_node = 'Root'):
        new_clade = Phylo.BaseTree.Clade(name = add_node, clades = [self.phylo_tree.root])
        self.phylo_tree.root = new_clade

    def root_with_duplication(self, node_name): # This function is restricted to change root location
        assert(self.is_duplication_node(node_name))
        duplication_clade = self.find_clade(node_name)
        assert(duplication_clade in self.phylo_tree.root.clades)
        if len(self.phylo_tree.root.clades) == 2:
            duplication_clade.clades.extend([clade for clade in self.phylo_tree.root.clades if clade.name != node_name])
        self.phylo_tree.root = duplication_clade

        # update json tree
        self.get_tree_json()

        
    def update_tree(self):
        for i in range(len(self.edge_list)):
            node1 = self.num_to_node[self.tree_json['row_nodes'][i]]
            node2 = self.num_to_node[self.tree_json['column_nodes'][i]]
            self.tree_json['edge_rate_scaling_factors'][i] = self.edge_to_blen[(node1, node2)]

    def unpack_x_rates(self, log_x_rates, Force_rates = None):  # TODO: Change it to fit general tree structure rather than cherry tree
        assert(len(log_x_rates) == len(self.edge_list))
        x_rates = np.exp(log_x_rates)

        if Force_rates != None:
            for i in Force_rates.keys():
                x_rates[i] = Force_rates[i]
        assert(len(x_rates) == len(self.edge_to_blen))

        for edge_it in range(len(self.edge_list)):
            self.edge_to_blen[self.edge_list[edge_it]] = x_rates[edge_it] 

        self.update_tree()


    def find_clade(self, clade_name):
        hit_clades = list(self.phylo_tree.find_clades(clade_name))
        assert(len(hit_clades) == 1)
        return hit_clades[0]

    def find_parent_clade(self, clade_name):
        child_clade = self.find_clade(clade_name)
        path = self.phylo_tree.get_path(child_clade)
        if len(path) > 1:
            father_clade = path[-2]
        else:
            father_clade = self.phylo_tree.root

        return father_clade

       
    def init_root_conf(self):
        self.node_to_conf[self.phylo_tree.root.name] = [[0, 1]]
        self.n_orlg = 1
    
    def get_configurations_for_path(self, terminal_node_name, node_to_pos):
        terminal_clade = self.find_clade(terminal_node_name)
        path = self.phylo_tree.get_path(terminal_clade)
        for clade in path:
            if clade.name in self.node_to_conf:
                continue
            elif self.is_duplication_node(clade.name):
                self.get_configuration_for_duplication_node(clade.name, node_to_pos[clade.name])
                self.visited_DL_nodes.append(clade.name)
            elif self.is_deletion_node(clade.name):
                self.get_configuration_for_deletion_node(clade.name, node_to_pos[clade.name])
                self.visited_DL_nodes.append(clade.name)
            elif self.is_speciation_node(clade.name):
                self.get_configuration_for_speciation_node(clade.name)
            elif self.is_terminal_node(clade.name):
                self.get_configuration_for_terminal_node(clade.name)
            else:
                print 'The node cannot be recognised!'
        

    def get_configurations(self):
        all_DP_nodes = [node for node in self.node_to_num if self.is_duplication_node(node) or self.is_deletion_node(node)]
        assert(all([node in self.node_to_pos for node in all_DP_nodes]))
        if not self.phylo_tree.root.name in self.node_to_conf:
            self.init_root_conf()
        for node in self.terminal_node_list:
            self.get_configurations_for_path(node, self.node_to_pos)

        assert(all([node in self.visited_DL_nodes for node in self.node_to_pos]))
        assert(self.is_configurations_same_size())
        self.n_js = len(self.node_to_conf[self.node_to_conf.keys()[0]])

    def is_duplication_node(self, node_name):
        return node_name[0] == 'D' and str.isdigit(node_name[1:])

    def is_deletion_node(self, node_name):
        return node_name[0] == 'L' and str.isdigit(node_name[1:])

    def is_speciation_node(self, node_name):
        return node_name[0] == 'N' and str.isdigit(node_name[1:])

    def is_terminal_node(self, node_name):
        return node_name in [clade.name for clade in self.phylo_tree.get_terminals()]
        
    def get_configuration_for_terminal_node(self, node_name):
        # A terminal node copies its configuration from its parent node
        self.copy_configuration_from_parent(node_name)

    def copy_configuration_from_parent(self, node_name):
        parent_clade = self.find_parent_clade(node_name)
        self.node_to_conf[node_name] = deepcopy(self.node_to_conf[parent_clade.name])

    def get_configuration_for_speciation_node(self, node_name):
        # A speciation node copies its configuration from its parent node
        self.copy_configuration_from_parent(node_name)

    def is_configurations_same_size(self):
        return len(set([len(self.node_to_conf[node]) for node in self.node_to_conf])) == 1

    def get_configuration_for_duplication_node(self, node_name, orlg_pos): # this is simplified version
        # There are cases this function cannot handle
        # For example, this duplication tree: ((3, 4)D2, (1, 2)D1)N0
        assert(not node_name in self.node_to_conf)
        parent_clade = self.find_parent_clade(node_name)
        assert(parent_clade.name in self.node_to_conf)
        old_configuration = self.node_to_conf[parent_clade.name]
        ortho_group_to_pos = divide_configuration(old_configuration)
        
        if type(orlg_pos) == int:
            old_orlg = ortho_group_to_pos['loc'][orlg_pos]
            
            assert(self.is_configurations_same_size())
            if len(ortho_group_to_pos['extent'][old_orlg]) == 1:# now only deal with case that the ortholog group occupies only one position
                pos = ortho_group_to_pos['extent'][old_orlg][0]
                assert(self.node_to_conf[parent_clade.name][pos][1]) # only extent lineage can give birth

                # Step 1, replace old ortholog group with two new groups
                new_orlg_1 = self.n_orlg
                new_orlg_2 = self.n_orlg + 1
                self.n_orlg += 2
                self.dup_events[old_orlg] = [new_orlg_1, new_orlg_2]
                self.node_to_dup[node_name] = [new_orlg_1, new_orlg_2]

                # Step 2, update all other configurations
                # They should be of same size as old_configuration
                for node in self.node_to_conf:
                    self.node_to_conf[node].insert(pos, deepcopy(self.node_to_conf[node][pos])) # duplicate the parent position
                
                # insert current node's configuration
                new_configuration = deepcopy(old_configuration)
                new_configuration[pos][0] = new_orlg_1
                new_configuration[pos + 1][0] = new_orlg_2
                self.node_to_conf[node_name] = new_configuration
            else:
                divided_positions = self.divide_positions(ortho_group_to_pos['extent'][old_orlg])
                assert(len(divided_positions) == 1) # Now only consider one continuous representation
                # TODO:implement more general case
                pos_list = divided_positions[0]

                # Step 1, replace old ortholog group with two new groups
                new_orlg_1 = self.n_orlg
                new_orlg_2 = self.n_orlg + 1
                self.n_orlg += 2
                self.dup_events[old_orlg] = [new_orlg_1, new_orlg_2]
                self.node_to_dup[node_name] = [new_orlg_1, new_orlg_2]

                # No need for Step 2: update all other configurations
                # They should stay of the same size as old_configuration

                # insert current node's configuration
                new_configuration = deepcopy(old_configuration)
                for i in range(len(pos_list)):
                    pos = pos_list[i]
                    if i < len(pos_list) / 2:
                        new_configuration[pos_list[i]][0] = new_orlg_1
                    else:
                        new_configuration[pos_list[i]][0] = new_orlg_2
                self.node_to_conf[node_name] = new_configuration
        elif type(orlg_pos) == list:
            old_orlg_list = [ortho_group_to_pos['loc'][i] for i in orlg_pos]
            sys.exit('TODO: Implement duplication events affecting multiple contiguous paralogs')
        else:
            sys.exit('Please check get_configuration_for_duplication_node() function in Tree class')

        
    def divide_positions(self, pos_list):
        results = []
        for k, g in groupby(enumerate(pos_list), lambda (i, x):i-x):
            results.append(map(itemgetter(1), g))
        return results

    def get_configuration_for_deletion_node(self, node_name, orlg_pos):
        parent_clade = self.find_parent_clade(node_name)
        assert(parent_clade.name in self.node_to_conf)
        new_configuration = deepcopy(self.node_to_conf[parent_clade.name])        
        ortho_group_to_pos = divide_configuration(new_configuration)
        deleted_orlg = ortho_group_to_pos['loc'][orlg_pos]
        for pos in ortho_group_to_pos['extent'][deleted_orlg]:
            assert(new_configuration[pos][1])  # the paralog should be alive before deletion
            new_configuration[pos][1] = 0
        self.node_to_conf[node_name] = new_configuration

    def divide_configuration(self, configuration):
        ortho_group_to_pos = dict(extent = {}, distinct = [], loc = [])
        # extent positions that represent same paralog (by the same ortholog group number) have to be in the same state
        # distinct positions don't change states, thus only track positions
        for pos in range(len(configuration)):
            if configuration[pos][1] == 1: # extent
                ortho_group = configuration[pos][0]
                if ortho_group in ortho_group_to_pos['extent']:
                    ortho_group_to_pos['extent'][ortho_group].append(pos)
                else:
                    ortho_group_to_pos['extent'][ortho_group] = [pos]
                    ortho_group_to_pos['loc'].append(ortho_group)
            elif configuration[pos][1] == 0: # distinct
                ortho_group_to_pos['distinct'].append(pos)

        return ortho_group_to_pos

    def __str__(self):  # overide for print function
        Phylo.draw_ascii(self.phylo_tree)
        for node in sorted(self.node_to_conf.keys()):
            print node, self.node_to_conf[node]
        print
        return 'Tree newick file: ' + self.newicktree + '\n' + \
               'Tree duplos file: ' + self.duploslist + '\n'
               


    
if __name__ == '__main__':
    wd = '/Users/xji3/GitFolders/xji3ST790/CourseProject/'
    tree_newick = wd + 'sim_tree.newick'
    DupLosList = wd + 'Sim_DupLost.txt'
    terminal_node_list = ['Out', 'A', 'B']
    node_to_pos = node_to_pos = {'D1':0, 'D2':0}
    tree = Tree(tree_newick, DupLosList, terminal_node_list, node_to_pos)
    self = tree
    tree.get_tree()

##    tree_newick = '../test/PrimateTest.newick'
##    DupLosList = '../test/PrimateTestDupLost.txt'
##    tree = Phylo.read( tree_newick, "newick")
##    terminal_node_list = ['Chinese_Tree_Shrew', 'Macaque', 'Olive_Baboon', 'Orangutan', 'Gorilla', 'Human']
##    test = Tree(tree_newick, DupLosList)
##    Phylo.draw_ascii(test.phylo_tree)
##    self = test
####    father_node_name = 'N1'
####    child_node_name = 'N3'
####    test.unpack_x_rates(np.log([0.1] * len(test.edge_list)))
####    print test.edge_to_blen
####    print
####
####    test.init_root_conf()
####    terminal_node_name = 'Chinese_Tree_Shrew'
####    terminal_clade = self.find_clade(terminal_node_name)
####    path = self.phylo_tree.get_path(terminal_clade)
####    print path, test.is_duplication_node(path[0].name)
####
####    print test.node_to_conf
####
######    test.get_configuration_for_duplication_node('D0', 0)
######    print test.node_to_conf, test.dup_events
######
######    test.get_configuration_for_duplication_node('D5', 1)
######    print test.node_to_conf, test.dup_events
######
######    test.get_configuration_for_deletion_node('L0', 1)
######    print test.node_to_conf, test.dup_events
####
##    node_to_pos = {'D1':0, 'D2':0, 'D3':1, 'D4':2, 'L1':2}
####
####    test.get_configurations_for_path('Chinese_Tree_Shrew', node_to_pos)
####    print test.node_to_conf, test.dup_events
####
####    test.get_configurations_for_path('Macaque', node_to_pos)
####    print test.node_to_conf, test.dup_events
##
##    test.get_configurations(terminal_node_list, node_to_pos)
##    print test.dup_events
##    for i in test.node_to_conf:
##        if i in terminal_node_list:
##            print i, test.node_to_conf[i]
##
##
##    tree_newick = '../test/YeastTree.newick'
##    DupLosList = '../test/YeastTestDupLost.txt'
##    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
##    test = Tree(tree_newick, DupLosList)
##    Phylo.draw_ascii(test.phylo_tree)
##    node_to_pos = {'D1':0}
##    test.get_configurations(terminal_node_list, node_to_pos)
##
##    for i in test.node_to_conf:
##        if i in terminal_node_list:
##            print i, test.node_to_conf[i]


##    tree_newick = '../test/Trigeneconv_ADH1Class_tree.newick'
##    DupLosList = '../test/Trigeneconv_ADH_DupLost.txt'
##    terminal_node_list = ['Baboon', 'Orangutan', 'Gorilla', 'Bonobo', 'Chimpanzee', 'Human']
##    node_to_pos = {'D1':0, 'D2':0}
##    
##    test = Tree(tree_newick, DupLosList, terminal_node_list, node_to_pos)
##    
##    Phylo.draw_ascii(test.phylo_tree)
##    self = test
##    for node in test.node_to_conf:
##        print node, test.node_to_conf[node]

##    tree_newick = '../test/PrimateTree.newick'
##    DupLosList = '../test/PrimateFullDupLost_P2.txt'
##    terminal_node_list = ['Chinese_Tree_Shrew', 'Bushbaby', 'Mouse_Lemur',
##                          'Tarsier', 'Marmoset', 'Vervet-AGM',
##                          'Olive_Baboon', 'Macaque', 'Gibbon',
##                          'Orangutan', 'Gorilla', 'Human']
##    node_to_pos = {'D1':0, 'D2':0, 'D3':0, 'D4':0, 'D5':0}
##    test = Tree(tree_newick, DupLosList, terminal_node_list, node_to_pos)
##
##    Phylo.draw_ascii(test.phylo_tree)
##    self = test
##    test.get_configurations()
##    for i in test.node_to_conf:
##        print i, test.node_to_conf[i]

#    
            
##    test.root_with_duplication('D1')
##    Phylo.draw_ascii(test.phylo_tree)
##    for node in test.node_to_conf:
##        print node, test.node_to_conf[node]

