# Xiang decided to seperate functions and the recording class so that he can easily vectorize his code
# xji3@ncsu.edu

# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Re-written of my previous CondonGeneconv class

from operator import mul
from itertools import product
from functools import partial
from copy import deepcopy
import os, sys

import numpy as np
#import networkx as nx
import scipy
import scipy.optimize
import scipy.sparse
import scipy.sparse.linalg

from Bio import Phylo
from Bio import SeqIO

#import cProfile
import jsonctmctree.ll, jsonctmctree.interface

##def get_HKYBasicRate(na, nb, pi, kappa):
##    if isTransition(na, nb):
##        return pi['ACGT'.index(nb)] * kappa
##    else:
##        return pi['ACGT'.index(nb)]

def get_HKYGeneconvRate(pair_from, pair_to, Qbasic, tau):
    na, nb = pair_from
    nc, nd = pair_to
    if (na != nc and nb!= nd) or (na == nc and nb == nd):
        return 0.0
    if na ==nc and nb != nd:
        Qb = Qbasic['ACGT'.index(nb), 'ACGT'.index(nd)]
        if na == nd:
            return Qb + tau
        else:
            return Qb
    if nb == nd and na != nc:
        Qb = Qbasic['ACGT'.index(na), 'ACGT'.index(nc)]
        if nb == nc:
            return Qb + tau
        else:
            return Qb
    print ('Warning: Check get_HKYGeneconvRate Func. You should not see this.')

def get_MG94BasicRate(ca, cb, pi, kappa, omega, codon_table):
    dif = [ii for ii in range(3) if ca[ii] != cb[ii]]
    ndiff = len(dif)
    if ndiff > 1:
        return 0
    elif ndiff == 0:
        print ('Please check your codon tables and make sure no redundancy')
        print (ca, cb)
        return 0
    else:
        na = ca[dif[0]]
        nb = cb[dif[0]]
        QbasicRate = pi['ACGT'.index(nb)]

        if isTransition(na, nb):
            QbasicRate *= kappa

        if isNonsynonymous(ca, cb, codon_table):
            QbasicRate *= omega

        return QbasicRate

def isTransition(na, nb):
    return (set([na, nb]) == set(['A', 'G']) or set([na, nb]) == set(['C', 'T']))

def isNonsynonymous(ca, cb, codon_table):
    return (codon_table[ca] != codon_table[cb])

#vec_get_MG94BasicRate = np.vectorize(get_MG94BasicRate, doc='Vectorized `get_MG94BasicRate`', excluded = ['pi', 'kappa', 'omega', 'codon_table'])

def get_MG94GeneconvRate(pair_from, pair_to, Qbasic, tau, codon_to_state):
    # pair_from = a string of length 6
    ca, cb = pair_from[:3], pair_from[3:]
    cc, cd = pair_to[:3], pair_to[3:]
    row = (codon_to_state[ca], codon_to_state[cb])
    col = (codon_to_state[cc], codon_to_state[cd])
    if ca != cc and cb != cd:
        return None
    if ca == cc and cb == cd: # do not deal with diagonal entries here
        return None

    if ca == cc:
        # cb != cd
        BasicRate = Qbasic[codon_to_state[cb], codon_to_state[cd]]
        if cd == ca:
            if isNonsynonymous(cb, cd, codon_table):
                additional_source = tau * omega
            else:
                additional_source = tau
            return [row, col, BasicRate + additional_source]
        else:
            return [row, col, BasicRate]
    else: # cb == cd
        # ca != cc
        BasicRate = Qbasic[codon_to_state[ca], codon_to_state[cc]]
        if cc == cb:
            if isNonsynonymous(ca, cc, codon_table):
                additional_source = tau * omega
            else:
                additional_source = tau
            return [row, col, BasicRate + additional_source]
        else:
            return [row, col, BasicRate]
    #print ('You should not see this')

## vec_get_MG94GeneconvRate = np.vectorize(get_MG94GeneconvRate, doc = 'Vectorized `get_MG94GeneconvRate`', excluded = ['pi', 'kappa', 'omega', 'codon_table', 'tau', 'codon_to_state'])

def get_x_clock_guess(edge_to_blen):
    # TODO: modify this to work on general tree topology
    leaf_branch = [edge for edge in edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
    out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
    internal_branch = [x for x in edge_to_blen.keys() if not x in leaf_branch]


    leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
    internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order

    Lr_reverse = []
    for i in range(len(internal_branch) - 1, 0, -1):
        Lr_reverse.append(edge_to_blen[leaf_branch[i]]/(edge_to_blen[leaf_branch[i]] + edge_to_blen[internal_branch[i]]))
    r0 = 2 / (edge_to_blen[out_group_branch] * (1 - Lr_reverse[-1]) / edge_to_blen[internal_branch[0]] + 1)
    Lr_reverse.append(r0)
    Lr_reverse.append(edge_to_blen[out_group_branch] / (2 - r0))
    return Lr_reverse

def read_newick(newick_file, post_dup = 'N1'):
    assert(os.path.isfile(newick_file))  # check if file exists
    tree = Phylo.read(newick_file, 'newick', rooted=True)

    # locate 1st post-duplication node
    post_dup_clade = tree.find_clades(post_dup).next()

    # if the root node is the 1st post-duplication node,
    # then add a duplication node before the root
    # inspired by root_with_outgroup() function in Phylo
    # http://biopython.org/DIST/docs/api/Bio.Phylo.BaseTree-pysrc.html
    if tree.root == post_dup_clade:
        new_root = tree.root.__class__(branch_length = tree.root.branch_length, name = 'dup')
        new_root.clades.insert(0, tree.root)
        tree.root = new_root
    
    ## from http://biopython.org/wiki/Phylo_cookbook
    allclades = list(tree.find_clades(order = 'level'))
    node_to_num = {n.name:i for i, n in enumerate(allclades)}
    
    edge_list = []
    for parent in tree.find_clades(terminal = False, order = 'level'):
        for child in parent.clades:
            edge_list.append((parent, child))

    tree_row = [node_to_num[na.name] for na, nb in edge_list]
    tree_col = [node_to_num[nb.name] for na, nb in edge_list]
    tree_process = [1 if post_dup_clade.is_parent_of(edge[1]) else 0 for edge in edge_list]

    out_tree = dict(
            row = tree_row,
            col = tree_col,
            process = tree_process,
            rate = np.ones(len(tree_row))
            )

    edge_list = [(clade_a.name, clade_b.name) for clade_a, clade_b in edge_list]
    return out_tree, edge_list, node_to_num
      

if __name__ == '__main__':
    newick_file = '/Users/xji3/GitFolders/Genconv/IGCexpansion/YeastTree_test.newick'

    post_dup = 'N1'
    tree = Phylo.read(newick_file, 'newick')

    # locate 1st post-duplication node
    post_dup_clade = tree.find_clades(post_dup).next()
    
    ## from http://biopython.org/wiki/Phylo_cookbook
    allclades = list(tree.find_clades(order = 'level'))
    node_to_num = {n.name:i for i, n in enumerate(allclades)}
    num_to_node = {node_to_num[n]:n for n in node_to_num.keys()}

    edge_list = []
    for parent in tree.find_clades(terminal = False, order = 'level'):
        for child in parent.clades:
            edge_list.append((parent, child))

    tree_row = [node_to_num[na.name] for na, nb in edge_list]
    tree_col = [node_to_num[nb.name] for na, nb in edge_list]
    tree_process = [1 if post_dup_clade.is_parent_of(edge[1]) else 0 for edge in edge_list]

    out_tree = dict(
            row = tree_row,
            col = tree_col,
            process = tree_process,
            rate = np.ones(len(tree_row))
            )                                       
    
