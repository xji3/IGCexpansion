# A separate class to represent tree structure
# Xiang Ji
# xji3@ncsu.edu
from Bio import Phylo
import networkx as nx
import os
import numpy as np

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

    def add_duplos_nodes(self):
        assert(os.path.isfile(self.duploslist))
        with open(self.duploslist, 'rb') as f:
            for line in f:
                items = line.split()
                if items:
                    branch = items[0]
                    father_node_name, child_node_name = branch.split('_')
                    for add_node in items[1:]:
                        self.add_node(father_node_name, child_node_name, add_node)
                        father_node_name = add_node
                    
    def add_node(self, father_node_name, child_node_name, add_node):
        hit_clades = list(self.phylo_tree.find_clades(child_node_name))
        assert(len(hit_clades) == 1)
        child_clade = hit_clades[0]

        father_clade = self.phylo_tree.get_path(child_clade)[-2]
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

if __name__ == '__main__':
    tree_newick = '../test/PrimateTest.newick'
    DupLosList = '../test/PrimateTestDupLost.txt'
    tree = Phylo.read( tree_newick, "newick")
    test = Tree(tree_newick, DupLosList)
    self = test
    father_node_name = 'N1'
    child_node_name = 'N3'
    Phylo.draw_ascii(self.phylo_tree)
    test.unpack_x_rates(np.log([0.1] * len(test.edge_list)))
    print test.edge_to_blen
    print
