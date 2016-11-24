# A separate class to represent tree structure
# Xiang Ji
# xji3@ncsu.edu
from Bio import Phylo
import networkx as nx
import os
import numpy as np
from copy import deepcopy

class Tree:
    def __init__(self, tree_newick, DupLosList):
        self.newicktree   = tree_newick   # newick tree file location
        self.duploslist   = DupLosList    # duplication loss nodes file location
        # Tree topology related variable
        # The tree for this project is fixed, but keep the read-in feature anyway
        self.phylo_tree   = None        # Phylo tree structure
        self.tree_json    = None        # store the tree dictionary used for json likelihood package parsing
        self.edge_to_blen = None        # dictionary store the unpacked tree branch length information {(node_from, node_to):blen}
        self.edge_list    = None        # kept all edges in the same order with x_rates
        self.node_to_num  = None        # dictionary used for translating tree info from self.edge_to_blen to self.tree
        self.num_to_node  = None        # dictionary used for translating tree info from self.tree to self.edge_to_blen

        self.node_to_conf = dict()      # A dictionary store configurations on each node
        # Speciation node starts with N, Duplication node with D, Loss node with L
        self.dup_events   = dict()      # A dictionary stores duplication events: ortholog group in key gives birth to the two ortholog groups in the content list

        self.n_orlg       = 0           # number of ortholog groups
        self.get_tree()



    def get_tree(self):
        tree = Phylo.read( self.newicktree, "newick")
        self.phylo_tree = tree.as_phyloxml(rooted = 'True')
        self.add_duplos_nodes()
        #set node number for nonterminal nodes and specify root node
        
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
        edge_list = sorted(edge_to_blen.keys(), key = lambda edge: edge[0])
         
        # Now setup self.tree dictionary
        tree_row = [self.node_to_num[na] for na, nb in edge_list]
        tree_col = [self.node_to_num[nb] for na, nb in edge_list]
        tree_process = [1 if e[0] == 'N0' and e[1] == 'N1' else 2 for e in edge_list]
        self.edge_list = edge_list

        self.tree_json = dict(
            row = tree_row,
            col = tree_col,
            process = tree_process,
            rate = np.ones(len(tree_row))
            )

        # TODO: need to change process part

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

    def update_tree(self):
        for i in range(len(self.tree_json['rate'])):
            node1 = self.num_to_node[self.tree_json['row'][i]]
            node2 = self.num_to_node[self.tree_json['col'][i]]
            self.tree_json['rate'][i] = self.edge_to_blen[(node1, node2)]

    def unpack_x_rates(self, log_x_rates, Force_rates = None):  # TODO: Change it to fit general tree structure rather than cherry tree
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
            elif self.is_deletion_node(clade.name):
                self.get_configuration_for_deletion_node(clade.name, node_to_pos[clade.name])
            elif self.is_speciation_node(clade.name):
                self.get_configuration_for_speciation_node(clade.name)
            elif self.is_terminal_node(clade.name):
                self.get_configuration_for_terminal_node(clade.name)
            else:
                print 'The node cannot be recognised!'

    def get_configurations(self, terminal_node_list, node_to_pos):
        self.init_root_conf()
        for node in terminal_node_list:
            self.get_configurations_for_path(node, node_to_pos)

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
        parent_clade = self.find_parent_clade(node_name)
        assert(parent_clade.name in self.node_to_conf)
        old_configuration = self.node_to_conf[parent_clade.name]
        ortho_group_to_pos = self.divide_configuration(old_configuration)
        old_orlg = sorted(ortho_group_to_pos['extent'].keys())[orlg_pos]
        
        assert(self.is_configurations_same_size())
        assert(len(ortho_group_to_pos['extent'][old_orlg]) == 1) # now only deal with case that the ortholog group occupies only one position
        pos = ortho_group_to_pos['extent'][old_orlg][0]
        assert(self.node_to_conf[parent_clade.name][pos][1]) # only extent lineage can give birth

        # Step 1, replace old ortholog group with two new groups
        new_orlg_1 = self.n_orlg
        new_orlg_2 = self.n_orlg + 1
        self.n_orlg += 2
        self.dup_events[old_orlg] = [new_orlg_1, new_orlg_2]

        # Step 2, update all other configurations
        # They should be of same size as old_configuration
        for node in self.node_to_conf:
            self.node_to_conf[node].insert(pos, deepcopy(self.node_to_conf[node][pos])) # duplicate the parent position
        
        # insert current node's configuration
        new_configuration = deepcopy(old_configuration)
        new_configuration[pos][0] = new_orlg_1
        new_configuration[pos + 1][0] = new_orlg_2
        self.node_to_conf[node_name] = new_configuration

    def get_configuration_for_deletion_node(self, node_name, orlg_pos):
        parent_clade = self.find_parent_clade(node_name)
        assert(parent_clade.name in self.node_to_conf)
        new_configuration = deepcopy(self.node_to_conf[parent_clade.name])        
        ortho_group_to_pos = self.divide_configuration(new_configuration)
        deleted_orlg = sorted(ortho_group_to_pos['extent'].keys())[orlg_pos]
        for pos in ortho_group_to_pos['extent'][deleted_orlg]:
            assert(new_configuration[pos][1])  # the paralog should be alive before deletion
            new_configuration[pos][1] = 0
        self.node_to_conf[node_name] = new_configuration

    def divide_configuration(self, configuration):
        ortho_group_to_pos = dict(extent = {}, distinct = [])
        # extent positions that represent same paralog (by the same ortholog group number) have to be in the same state
        # distinct positions don't change states, thus only track positions
        for pos in range(len(configuration)):
            if configuration[pos][1] == 1: # extent
                ortho_group = configuration[pos][0]
                if ortho_group in ortho_group_to_pos['extent']:
                    ortho_group_to_pos['extent'][ortho_group].append(pos)
                else:
                    ortho_group_to_pos['extent'][ortho_group] = [pos]
            elif configuration[pos][1] == 0: # distinct
                ortho_group_to_pos['distinct'].append(pos)

        return ortho_group_to_pos        


    
if __name__ == '__main__':
    tree_newick = '../test/PrimateTest.newick'
    DupLosList = '../test/PrimateTestDupLost.txt'
    tree = Phylo.read( tree_newick, "newick")
    test = Tree(tree_newick, DupLosList)
    Phylo.draw_ascii(test.phylo_tree)
    self = test
##    father_node_name = 'N1'
##    child_node_name = 'N3'
##    test.unpack_x_rates(np.log([0.1] * len(test.edge_list)))
##    print test.edge_to_blen
##    print
##
##    test.init_root_conf()
##    terminal_node_name = 'Chinese_Tree_Shrew'
##    terminal_clade = self.find_clade(terminal_node_name)
##    path = self.phylo_tree.get_path(terminal_clade)
##    print path, test.is_duplication_node(path[0].name)
##
##    print test.node_to_conf
##
####    test.get_configuration_for_duplication_node('D0', 0)
####    print test.node_to_conf, test.dup_events
####
####    test.get_configuration_for_duplication_node('D5', 1)
####    print test.node_to_conf, test.dup_events
####
####    test.get_configuration_for_deletion_node('L0', 1)
####    print test.node_to_conf, test.dup_events
##
    node_to_pos = {'D1':0, 'D2':0, 'D3':1, 'D4':2, 'L1':2}
##
##    test.get_configurations_for_path('Chinese_Tree_Shrew', node_to_pos)
##    print test.node_to_conf, test.dup_events
##
##    test.get_configurations_for_path('Macaque', node_to_pos)
##    print test.node_to_conf, test.dup_events

    terminal_node_list = ['Chinese_Tree_Shrew', 'Macaque', 'Olive_Baboon', 'Orangutan', 'Gorilla', 'Human']
    test.get_configurations(terminal_node_list, node_to_pos)
    print test.dup_events
    for i in test.node_to_conf:
        if i in terminal_node_list:
            print i, test.node_to_conf[i]

 
