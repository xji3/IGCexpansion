# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Re-write of my previous CodonGeneconv class
# commit number: Oct 22nd, 2014 for old package
# cb1ba60ee2b57d6703cd9a3987000c2fd4dd68a5
# commit number: Dec 17th, 2014 for new package
# 33e393a973161e3a29149e82bfda23882b5826f3
from __future__ import print_function, absolute_import
from IGCexpansion.CodonGeneconFunc import *
import argparse
#from jsonctmctree.extras import optimize_em
import ast

#import matplotlib.pyplot as plt

class ReCodonGeneconv:
    def __init__(self, tree_newick, alignment, paralog, Model = 'MG94', nnsites = None, clock = False, Force = None, save_path = './save/', save_name = None):
        self.newicktree  = tree_newick  # newick tree file loc
        self.seqloc      = alignment    # multiple sequence alignment, now need to remove gap before-hand
        self.paralog     = paralog      # parlaog list
        self.nsites      = nnsites      # number of sites in the alignment used for calculation
        self.Model       = Model
        self.ll          = 0.0          # Store current log-likelihood
        self.Force       = Force        # parameter constraints only works on self.x not on x_clock which should be translated into self.x first
        self.clock       = clock        # molecular clock control
        self.save_path   = save_path    # location for auto-save files
        self.save_name   = save_name    # save file name
        self.auto_save   = 0            # auto save control

        self.logzero     = -15.0        # used to avoid log(0), replace log(0) with -15
        self.infinity    = 1e6          # used to avoid -inf in gradiance calculation of the clock case
        self.minlogblen  = -9.0         # log value, used for bond constraint on branch length estimates in get_mle() function

        # Tree topology related variable
        self.tree         = None        # store the tree dictionary used for json likelihood package parsing
        self.edge_to_blen = None        # dictionary store the unpacked tree branch length information {(node_from, node_to):blen}
        self.edge_list    = None        # kept all edges in the same order with x_rates
        self.node_to_num  = None        # dictionary used for translating tree info from self.edge_to_blen to self.tree
        self.num_to_node  = None        # dictionary used for translating tree info from self.tree to self.edge_to_blen

        # Constants for Sequence operations
        bases = 'tcag'.upper()
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        
        self.nt_to_state    = {a:i for (i, a) in enumerate('ACGT')}
        self.codon_table    = dict(zip(codons, amino_acids))
        self.codon_nonstop  = [a for a in self.codon_table.keys() if not self.codon_table[a]=='*']
        self.codon_to_state = {a.upper() : i for (i, a) in enumerate(self.codon_nonstop)}
        self.pair_to_state  = {pair:i for i, pair in enumerate(product(self.codon_nonstop, repeat = 2))}

        # Tip data related variable
        self.name_to_seq      = None    # dictionary store sequences
        self.observable_names = None    # list of extent species + paralog name ( = self.name_to_seq.keys())
        self.observable_nodes = None    # list of extent species numbers (node_to_num)
        self.observable_axes  = None    # list of paralog numbers  
        self.iid_observations = None    # list of multivariate states


        # Rate matrix related variable
        self.x_process      = None      # values of process parameters (untransformed, log, or exp(-x))
        self.x_rates        = None      # values of blen (untransformed, log, or exp(-x))
        self.x              = None      # x_process + x_rates
        self.x_Lr           = None      # values of clock blen parameters
        self.x_clock        = None      # x_process + Lr
        self.pi             = None      # real values
        self.kappa          = 1.2       # real values
        self.omega          = 0.9       # real values
        self.tau            = 1.4       # real values

        self.processes      = None      # list of basic and geneconv rate matrices. Each matrix is a dictionary used for json parsing

        self.scene_ll       = None      # used for lnL calculation

        # Prior distribution on the root
        self.prior_feasible_states  = None
        self.prior_distribution     = None

        # Expected_Geneconv Events
        self.GeneconvTransRed  = None    # dictionary of Geneconv transition matrix used for json parsing
        self.ExpectedGeneconv  = None    # dictionary storing expected number of geneconv events on each branch
        self.ExpectedDwellTime = None    # dictionary storing expected total dwell time of heterogeneous states of each branch same order as self.edge_list

        # Initialize all parameters
        self.initialize_parameters()
        
        
    def initialize_parameters(self):
        self.get_tree()
        self.get_data()
        self.get_initial_x_process()
        if self.save_name == None:
            save_file = self.gen_save_file_name()
        else:
            save_file = self.save_name
        if os.path.isfile(save_file):  # if the save txt file exists and not empty, then read in parameter values
            if os.stat(save_file).st_size > 0:
                self.initialize_by_save(save_file)
                print ('Successfully loaded parameter value from ' + save_file)
                
    def get_tree(self):
        self.tree, self.edge_list, self.node_to_num = read_newick(self.newicktree)
        self.num_to_node = {self.node_to_num[i]:i for i in self.node_to_num}
        self.edge_to_blen = {edge:1.0 for edge in self.edge_list}

    def get_tree_old(self):
        tree = Phylo.read( self.newicktree, "newick")
        #set node number for nonterminal nodes and specify root node
        numInternalNode = 0
        for clade in tree.get_nonterminals():
            clade.name = 'N' + str(numInternalNode)
            numInternalNode += 1
        tree_phy = tree.as_phyloxml(rooted = 'True')
        tree_nx = Phylo.to_networkx(tree_phy)

        triples = ((u.name, v.name, d['weight']) for (u, v, d) in tree_nx.edges(data = True)) # data = True to have the blen as 'weight'
        T = nx.DiGraph()
        edge_to_blen = {}
        for va, vb, blen in triples:
            edge = (va, vb)
            T.add_edge(*edge)
            edge_to_blen[edge] = blen

        self.edge_to_blen = edge_to_blen

        # Now assign node_to_num
        leaves = set(v for v, degree in T.degree().items() if degree == 1)
        internal_nodes = set(list(T)).difference(leaves)
        node_names = list(internal_nodes) + list(leaves)
        self.node_to_num = {n:i for i, n in enumerate(node_names)}
        self.num_to_node = {self.node_to_num[i]:i for i in self.node_to_num}

        # Prepare for generating self.tree so that it has same order as the self.x_process
        nEdge = len(self.edge_to_blen)  # number of edges
        l = nEdge / 2 + 1               # number of leaves
        k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

        leaf_branch = [edge for edge in self.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
        out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
        internal_branch = [x for x in self.edge_to_blen.keys() if not x in leaf_branch]
        assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

        leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        edge_list = []
        for i in range(len(internal_branch)):
            edge_list.append(internal_branch[i])
            edge_list.append(leaf_branch[i])
        for j in range(len(leaf_branch[i + 1:])):
            edge_list.append(leaf_branch[i + 1 + j])
            
        # Now setup self.tree dictionary
        tree_row = [self.node_to_num[na] for na, nb in edge_list]
        tree_col = [self.node_to_num[nb] for na, nb in edge_list]
        tree_process = [0 if e[0] == 'N0' and e[1] != 'N1' else 1 for e in edge_list]
        self.edge_list = edge_list

        self.tree = dict(
            row = tree_row,
            col = tree_col,
            process = tree_process,
            rate = np.ones(len(tree_row))
            )

    def nts_to_codons(self):
        for name in self.name_to_seq.keys():
            assert(len(self.name_to_seq[name]) % 3 == 0)
            tmp_seq = [self.name_to_seq[name][3 * j : 3 * j + 3] for j in range(len(self.name_to_seq[name]) / 3 )]
            self.name_to_seq[name] = tmp_seq
       
    def get_data(self):
        seq_dict = SeqIO.to_dict(SeqIO.parse( self.seqloc, "fasta" ))
        self.name_to_seq = {name:str(seq_dict[name].seq) for name in seq_dict.keys()}
        
        if self.Model == 'MG94':
            # Convert from nucleotide sequences to codon sequences.
            self.nts_to_codons()
            obs_to_state = self.codon_to_state
        else:
            obs_to_state = self.nt_to_state        

        # change the number of sites for calculation if requested
        if self.nsites is None:
            self.nsites = len(self.name_to_seq[self.name_to_seq.keys()[0]])
        else:
            for name in self.name_to_seq:
                self.name_to_seq[name] = self.name_to_seq[name][: self.nsites]
        print ('number of sites to be analyzed: ', self.nsites)

        # assign observable parameters
        suffix_len = len(self.paralog[0])
        self.observable_names = [n for n in self.name_to_seq.keys() if n[:-suffix_len] in self.node_to_num.keys()]
        paralog_len = [len(a) for a in self.paralog]
        assert(paralog_len[1:] == paralog_len[:-1])  # check if all paralog names have same length
        suffix_to_axis = {n:i for (i, n) in enumerate(list(set(self.paralog))) }
        self.observable_nodes = [self.node_to_num[n[:-suffix_len]] for n in self.observable_names]
        self.observable_axes = [suffix_to_axis[s[-suffix_len:]] for s in self.observable_names]
        
        # Now convert alignment into state list
        iid_observations = []
        for site in range(self.nsites):
            observations = []
            for name in self.observable_names:
                observation = obs_to_state[self.name_to_seq[name][site]]
                observations.append(observation)
            iid_observations.append(observations)
        self.iid_observations = iid_observations

    def get_initial_x_process(self, transformation = 'log'):
        
        count = np.array([0, 0, 0, 0], dtype = float) # count for A, C, G, T in all seq
        for name in self.name_to_seq:
            for i in range(4):
                count[i] += ''.join(self.name_to_seq[name]).count('ACGT'[i])
        count = count / count.sum()

        if self.Model == 'MG94':
            # x_process[] = %AG, %A, %C, kappa, omega, tau
            self.x_process = np.log(np.array([count[0] + count[2], count[0] / (count[0] + count[2]), count[1] / (count[1] + count[3]),
                                  self.kappa, self.omega, self.tau]))
        elif self.Model == 'HKY':
            # x_process[] = %AG, %A, %C, kappa, tau
            self.omega = 1.0
            self.x_process = np.log(np.array([count[0] + count[2], count[0] / (count[0] + count[2]), count[1] / (count[1] + count[3]),
                                  self.kappa, self.tau]))

        self.x_rates = np.log(np.array([ 0.1 * self.edge_to_blen[edge] for edge in self.edge_to_blen.keys()]))

        if transformation == 'log':    
            self.x = np.concatenate((self.x_process, self.x_rates))
        elif transformation == 'None':
            self.x_process = np.exp(self.x_process)
            self.x_rates = np.exp(self.x_rates)
        elif transformation == 'Exp_Neg':
            self.x_process = np.exp(-np.exp(self.x_process))
            self.x_rates = np.exp(-np.exp(self.x_rates))
        self.x = np.concatenate((self.x_process, self.x_rates))

        if self.clock:   # set-up x_clock if it's a clock model
            l = len(self.edge_to_blen) / 2 + 1               # number of leaves
            self.x_Lr = np.log(np.ones((l)) * 0.6)

            if transformation == 'log':
                self.x_clock = np.concatenate((self.x_process, self.x_Lr))
            elif transformation == 'None':
                self.x_Lr = np.exp(self.x_Lr)
            elif transformation == 'Exp_Neg':
                self.x_Lr = np.exp(-np.exp(self.x_Lr))
            self.x_clock = np.concatenate((self.x_process, self.x_Lr))
            self.unpack_x_clock(transformation = transformation)

        self.update_by_x(transformation = transformation)
        
    def update_by_x_clock(self, x_clock = None, transformation = 'log'):
        if not x_clock == None:
            self.x_clock = x_clock
        self.unpack_x_clock(transformation = transformation)
        self.update_by_x(transformation = transformation)
        
    def unpack_x_clock(self, transformation):
        assert(self.clock)
        nEdge = len(self.edge_to_blen)  # number of edges
        l = nEdge / 2 + 1               # number of leaves
        k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.
        if transformation == 'log':
            self.x_process, self.x_Lr = self.x_clock[:-l], np.exp(self.x_clock[-l:])
        elif transformation == 'None':
            self.x_process, self.x_Lr = self.x_clock[:-l], self.x_clock[-l:]
        elif transformation == 'Exp_Neg':
            self.x_process, self.x_Lr = self.x_clock[:-l], self.x_clock[-l:]
            self.x_clock[0] = - np.log(self.x_clock[0])

        # Now update self.x by using self.x_clock
        leaf_branch = [edge for edge in self.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
        out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
        internal_branch = [x for x in self.edge_to_blen.keys() if not x in leaf_branch]
        assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

        leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order

        # Now update blen with fixed order:
        # Always start from the root and internal-tip branch first
        for i in range(len(internal_branch)):
            edge = internal_branch[i]
            self.x_rates[2 * i] = self.blen_from_clock(edge)
            edge = leaf_branch[i]
            self.x_rates[2 * i + 1] = self.blen_from_clock(edge)
        for j in range(len(leaf_branch[i + 1:])):
            edge = leaf_branch[i + 1 + j]
            self.x_rates[ - len(leaf_branch[i + 1:]) + j] = self.blen_from_clock(edge)
        # update self.x so that all parameters can be updated by update_by_x
        if transformation == 'log':
            self.x_rates = np.array([ np.log(rate) if rate > 0 else self.logzero for rate in self.x_rates])
            self.x = np.concatenate((self.x_process, self.x_rates))
        elif transformation == 'None':
            self.x = np.concatenate((self.x_process, self.x_rates))
        elif transformation == 'Exp_Neg':
            self.x = np.concatenate((self.x_process, np.exp(-self.x_rates)))
        
        
    def blen_from_clock(self, edge):
        assert(edge in self.edge_to_blen.keys())
        if edge[0] == 'N0':
            if str.isdigit(edge[1][1:]):  # (N0, N1) branch
                return self.x_Lr[0] * self.x_Lr[1] * (1 - self.x_Lr[2])
            else:
                return self.x_Lr[0] * (2 - self.x_Lr[1])

        else:
            tmp_k = int(edge[0][1:])
            if str.isdigit( edge[1][1:] ): # ( N_temp_k, N_temp_k+1 ) branch
                return reduce( mul, self.x_Lr[: (tmp_k + 2)], 1)  * (1 - self.x_Lr[tmp_k + 2])
            else:  # ( N_temp_k, leaf ) branch
                return reduce( mul, self.x_Lr[: (tmp_k + 2)], 1)
        
        
    def update_by_x(self, x = None, transformation = 'log'):
        k = len(self.edge_to_blen)
        if not x == None:
            self.x = x
        self.x_process, self.x_rates = self.x[:-k], self.x[-k:]
        Force_process = None
        Force_rates = None
        if self.Force != None:
            Force_process = {i:self.Force[i] for i in self.Force.keys() if i < len(self.x) - k}
            Force_rates = {(i-len(self.x_process)):self.Force[i] for i in self.Force.keys() if not i < len(self.x) - k}
        self.unpack_x_process(Force_process = Force_process, transformation = transformation)
        self.unpack_x_rates(Force_rates = Force_rates, transformation = transformation)

    def unpack_x_process(self, transformation, Force_process = None):
        if transformation == 'log':
            x_process = np.exp(self.x_process)
        elif transformation == 'None':
            x_process = self.x_process
        elif transformation == 'Exp_Neg':
            x_process = x_process = np.concatenate((self.x_process[:3], -np.log(self.x_process[3:])))
            

        if Force_process != None:
            for i in Force_process.keys():
                x_process[i] = Force_process[i]

        if self.Model == 'MG94':
            # x_process[] = %AG, %A, %C, kappa, tau, omega
            assert(len(self.x_process) == 6)
            
            pi_a = x_process[0] * x_process[1]
            pi_c = (1 - x_process[0]) * x_process[2]
            pi_g = x_process[0] * (1 - x_process[1])
            pi_t = (1 - x_process[0]) * (1 - x_process[2])
            self.pi = [pi_a, pi_c, pi_g, pi_t]
            self.kappa = x_process[3]
            self.omega = x_process[4]
            self.tau = x_process[5]
        elif self.Model == 'HKY':
            # x_process[] = %AG, %A, %C, kappa, tau
            assert(len(self.x_process) == 5)
            pi_a = x_process[0] * x_process[1]
            pi_c = (1 - x_process[0]) * x_process[2]
            pi_g = x_process[0] * (1 - x_process[1])
            pi_t = (1 - x_process[0]) * (1 - x_process[2])
            self.pi = [pi_a, pi_c, pi_g, pi_t]
            self.kappa = x_process[3]
            self.tau = x_process[4]

        # Now update the prior distribution
        self.get_prior()

        # Now update processes (Rate matrices)
        self.get_processes()        

    def get_prior(self):
        if self.Model == 'MG94':
            self.prior_feasible_states = [(self.codon_to_state[codon], self.codon_to_state[codon]) for codon in self.codon_nonstop]
            distn = [ reduce(mul, [self.pi['ACGT'.index(b)]  for b in codon], 1) for codon in self.codon_nonstop ]
            distn = np.array(distn) / sum(distn)
        elif self.Model == 'HKY':
            self.prior_feasible_states = [(self.nt_to_state[nt], self.nt_to_state[nt]) for nt in 'ACGT']
            distn = [ self.pi['ACGT'.index(nt)] for nt in 'ACGT' ]
            distn = np.array(distn) / sum(distn)
        self.prior_distribution = distn

    def get_processes(self):
        if self.Model == 'MG94':
            self.processes = self.get_MG94Geneconv_and_MG94()
        elif self.Model == 'HKY':
            self.processes = self.get_HKYGeneconv()

    def get_MG94Geneconv_and_MG94(self):
        Qbasic = self.get_MG94Basic()
        row = []
        col = []
        rate_geneconv = []
        rate_basic = []

        for i, pair in enumerate(product(self.codon_nonstop, repeat = 2)):
            # use ca, cb, cc to denote codon_a, codon_b, codon_c, where cc != ca, cc != cb
            ca, cb = pair
            sa = self.codon_to_state[ca]
            sb = self.codon_to_state[cb]
            if ca != cb:
                for cc in self.codon_nonstop:
                    if cc == ca or cc == cb:
                        continue
                    sc = self.codon_to_state[cc]
                    # (ca, cb) to (ca, cc)
                    Qb = Qbasic[sb, sc]
                    if Qb != 0:
                        row.append((sa, sb))
                        col.append((sa, sc))
                        rate_geneconv.append(Qb)
                        rate_basic.append(0.0)

                    # (ca, cb) to (cc, cb)
                    Qb = Qbasic[sa, sc]
                    if Qb != 0:
                        row.append((sa, sb))
                        col.append((sc, sb))
                        rate_geneconv.append(Qb)
                        rate_basic.append(0.0)

                        
                # (ca, cb) to (ca, ca)
                row.append((sa, sb))
                col.append((sa, sa))
                Qb = Qbasic[sb, sa]
                if isNonsynonymous(cb, ca, self.codon_table):
                    Tgeneconv = self.tau * self.omega
                else:
                    Tgeneconv = self.tau
                rate_geneconv.append(Qb + Tgeneconv)
                rate_basic.append(0.0)
                
                # (ca, cb) to (cb, cb)
                row.append((sa, sb))
                col.append((sb, sb))
                Qb = Qbasic[sa, sb]
                rate_geneconv.append(Qb + Tgeneconv)
                rate_basic.append(0.0)

            else:
                for cc in self.codon_nonstop:
                    if cc == ca:
                        continue
                    sc = self.codon_to_state[cc]

                    # (ca, ca) to (ca,  cc)
                    Qb = Qbasic[sa, sc]
                    if Qb != 0:
                        row.append((sa, sb))
                        col.append((sa, sc))
                        rate_geneconv.append(Qb)
                        rate_basic.append(0.0)
                    # (ca, ca) to (cc, ca)
                        row.append((sa, sb))
                        col.append((sc, sa))
                        rate_geneconv.append(Qb)
                        rate_basic.append(0.0)

                    # (ca, ca) to (cc, cc)
                        row.append((sa, sb))
                        col.append((sc, sc))
                        rate_geneconv.append(0.0)
                        rate_basic.append(Qb)
                
        process_geneconv = dict(
            row = row,
            col = col,
            rate = rate_geneconv
            )
        process_basic = dict(
            row = row,
            col = col,
            rate = rate_basic
            )
        return [process_basic, process_geneconv]

    def get_MG94Basic(self):
        Qbasic = np.zeros((61, 61), dtype = float)
        for ca in self.codon_nonstop:
            for cb in self.codon_nonstop:
                if ca == cb:
                    continue
                Qbasic[self.codon_to_state[ca], self.codon_to_state[cb]] = get_MG94BasicRate(ca, cb, self.pi, self.kappa, self.omega, self.codon_table)
        expected_rate = np.dot(self.prior_distribution, Qbasic.sum(axis = 1))
        Qbasic = Qbasic / expected_rate
        return Qbasic

    def get_HKYBasic(self):
        Qbasic = np.array([
            [0, 1.0, self.kappa, 1.0],
            [1.0, 0, 1.0, self.kappa],
            [self.kappa, 1.0, 0, 1.0],
            [1.0, self.kappa, 1.0, 0],
            ]) * np.array(self.pi)
        expected_rate = np.dot(self.prior_distribution, Qbasic.sum(axis = 1))
        Qbasic = Qbasic / expected_rate
        return Qbasic
        
    
    def get_HKYGeneconv(self):
        #print ('tau = ', self.tau)
        Qbasic = self.get_HKYBasic()
        row = []
        col = []
        rate_geneconv = []
        rate_basic = []

        for i, pair_from in enumerate(product('ACGT', repeat = 2)):
            na, nb = pair_from
            sa = self.nt_to_state[na]
            sb = self.nt_to_state[nb]
            for j, pair_to in enumerate(product('ACGT', repeat = 2)):
                nc, nd = pair_to
                sc = self.nt_to_state[nc]
                sd = self.nt_to_state[nd]
                if i == j:
                    continue
                GeneconvRate = get_HKYGeneconvRate(pair_from, pair_to, Qbasic, self.tau)
                if GeneconvRate != 0.0:
                    row.append((sa, sb))
                    col.append((sc, sd))
                    rate_geneconv.append(GeneconvRate)
                    rate_basic.append(0.0)
                if na == nb and nc == nd:
                    row.append((sa, sb))
                    col.append((sc, sd))
                    rate_geneconv.append(GeneconvRate)
                    rate_basic.append(Qbasic['ACGT'.index(na), 'ACGT'.index(nc)])

        process_geneconv = dict(
            row = row,
            col = col,
            rate = rate_geneconv
            )
        process_basic = dict(
            row = row,
            col = col,
            rate = rate_basic
            )
        # process_basic is for HKY_Basic which is equivalent to 4by4 rate matrix
        return [process_basic, process_geneconv]
    
    def unpack_x_rates(self, transformation, Force_rates = None):  # TODO: Change it to fit general tree structure rather than cherry tree
        if transformation == 'log':
            x_rates = np.exp(self.x_rates)
        elif transformation == 'None':
            x_rates = self.x_rates
        elif transformation == 'Exp_Neg':
            x_rates = -np.log(self.x_rates)
            

        if Force_rates != None:
            for i in Force_rates.keys():
                x_rates[i] = Force_rates[i]
        assert(len(x_rates) == len(self.edge_to_blen))

        for edge_it in range(len(self.edge_list)):
            self.edge_to_blen[self.edge_list[edge_it]] = x_rates[edge_it] 

        self.update_tree()

    def update_tree(self):
        for i in range(len(self.tree['rate'])):
            node1 = self.num_to_node[self.tree['row'][i]]
            node2 = self.num_to_node[self.tree['col'][i]]
            self.tree['rate'][i] = self.edge_to_blen[(node1, node2)]

    def _loglikelihood(self, store = True, edge_derivative = False):
        '''
        Modified from Alex's objective_and_gradient function in ctmcaas/adv-log-likelihoods/mle_geneconv_common.py
        '''
        if self.Model == 'MG94':
            state_space_shape = [61, 61]
        elif self.Model == 'HKY':
            state_space_shape = [4, 4]

        # prepare some extra parameters for the json interface
        if edge_derivative:
            requested_derivatives = list(range(k))
        else:
            requested_derivatives = []
            
        site_weights = np.ones(self.nsites)

        # prepare the input for the json interface
        data = dict(
            site_weights = site_weights,
            requested_derivatives = requested_derivatives,
            node_count = len(self.edge_to_blen) + 1,
            state_space_shape = state_space_shape,
            process_count = len(self.processes),
            processes = self.processes,
            tree = self.tree,
            prior_feasible_states = self.prior_feasible_states,
            prior_distribution = self.prior_distribution,
            observable_nodes = self.observable_nodes,
            observable_axes = self.observable_axes,
            iid_observations = self.iid_observations
            )
        j_ll = jsonctmctree.ll.process_json_in(data)

        status = j_ll['status']
        feasibility = j_ll['feasibility']

        if status != 'success' or not feasibility:
            print ('results:')
            print (j_ll)
            print ()
            raise Exception('Encountered some problem in the calculation of log likelihood and its derivatives')

        ll, edge_derivs = j_ll['log_likelihood'], j_ll['edge_derivatives']
        self.ll = ll

        return ll, edge_derivs


    def _loglikelihood2(self, store = True, edge_derivative = False):
        '''
        Modified from Alex's objective_and_gradient function in ctmcaas/adv-log-likelihoods/mle_geneconv_common.py
        '''
        if store:
            self.scene_ll = self.get_scene()
            scene = self.scene_ll
        else:
            scene = self.get_scene()
        
        log_likelihood_request = {'property':'snnlogl'}
        derivatives_request = {'property':'sdnderi'}
        if edge_derivative:
            requests = [log_likelihood_request, derivatives_request]
        else:
            requests = [log_likelihood_request]
        j_in = {
            'scene' : self.scene_ll,
            'requests' : requests
            }
        j_out = jsonctmctree.interface.process_json_in(j_in)

        status = j_out['status']
    
        ll = j_out['responses'][0]
        self.ll = ll
        if edge_derivative:
            edge_derivs = j_out['responses'][1]
        else:
            edge_derivs = []

        return ll, edge_derivs
    
    def get_scene(self):
        if self.Model == 'MG94':
            state_space_shape = [61, 61]
        elif self.Model == 'HKY':
            state_space_shape = [4, 4]
        process_definitions = [{'row_states':i['row'], 'column_states':i['col'], 'transition_rates':i['rate']} for i in self.processes]
        scene = dict(
            node_count = len(self.edge_to_blen) + 1,
            process_count = len(self.processes),
            state_space_shape = state_space_shape,
            tree = {
                'row_nodes' : self.tree['row'],
                'column_nodes' : self.tree['col'],
                'edge_rate_scaling_factors' : self.tree['rate'],
                'edge_processes' : self.tree['process']
                },
            root_prior = {'states':self.prior_feasible_states,
                          'probabilities': self.prior_distribution},
            process_definitions = process_definitions,
            observed_data = {
                'nodes':self.observable_nodes,
                'variables':self.observable_axes,
                'iid_observations':self.iid_observations
                }            
            )
        return scene

    def loglikelihood_and_gradient(self, package = 'new', display = False):
        '''
        Modified from Alex's objective_and_gradient function in ctmcaas/adv-log-likelihoods/mle_geneconv_common.py
        '''
        self.update_by_x()
        delta = 1e-8
        x = deepcopy(self.x)  # store the current x array
        if package == 'new':
            fn = self._loglikelihood2
        else:
            fn = self._loglikelihood

        ll, edge_derivs = fn(edge_derivative = True)
        
        m = len(self.x) - len(self.edge_to_blen)

        # use finite differences to estimate derivatives with respect to these parameters
        other_derivs = []
        
        for i in range(m):
            if self.Force != None:
                if i in self.Force.keys():  # check here
                    other_derivs.append(0.0)
                    continue
            x_plus_delta = np.array(self.x)
            x_plus_delta[i] += delta
            self.update_by_x(x_plus_delta)
            ll_delta, _ = fn(store = True, edge_derivative = False)
            d_estimate = (ll_delta - ll) / delta           
            other_derivs.append(d_estimate)
            # restore self.x
            self.update_by_x(x)
        other_derivs = np.array(other_derivs)
        if display:
            print ('log likelihood = ', ll)
            print ('Edge derivatives = ', edge_derivs)
            print ('other derivatives:', other_derivs)
            print ('Current x array = ', self.x)

        self.ll = ll
        f = -ll
        g = -np.concatenate((other_derivs, edge_derivs))
        return f, g

    def loglikelihood_and_gradient2(self, package = 'new', display = False):
        '''
        Modified from Alex's objective_and_gradient function in ctmcaas/adv-log-likelihoods/mle_geneconv_common.py
        '''
        self.update_by_x()
        delta = 1e-8
        x = deepcopy(self.x)  # store the current x array
        if package == 'new':
            fn = self._loglikelihood2
        else:
            fn = self._loglikelihood

        ll, edge_derivs = fn(edge_derivative = True)
        
        m = len(self.x) - len(self.edge_to_blen)

        # use finite differences to estimate derivatives with respect to these parameters
        other_derivs = []
        
        for i in range(m):
            if self.Force != None:
                if i in self.Force.keys():  # check here
                    other_derivs.append(0.0)
                    continue
            x_plus_delta = np.array(self.x)
            x_plus_delta[i] += delta / 2.0
            self.update_by_x(x_plus_delta)
            ll_delta_plus, _ = fn(store = True, edge_derivative = False)
            x_plus_delta[i] -= delta
            self.update_by_x(x_plus_delta)
            ll_delta_minus, _ = fn(store = True, edge_derivative = False)
            x_plus_delta[i] += delta / 2.0
            d_estimate = (ll_delta_plus - ll_delta_minus) / delta           
            other_derivs.append(d_estimate)
            # restore self.x
            self.update_by_x(x)
        other_derivs = np.array(other_derivs)
        if display:
            print ('log likelihood = ', ll)
            print ('Edge derivatives = ', edge_derivs)
            print ('other derivatives:', other_derivs)
            print ('Current x array = ', self.x)

        self.ll = ll
        f = -ll
        g = -np.concatenate((other_derivs, edge_derivs))
        return f, g

    def objective_and_gradient(self, display, x):
        self.update_by_x(x)
        f, g = self.loglikelihood_and_gradient(display = display)
        self.auto_save += 1
        if self.auto_save == 5:
            self.save_x()
            self.auto_save = 0
        return f, g

    def Clock_wrap(self, display, x_clock):
        assert(self.clock)
        self.update_by_x_clock(x_clock)

        f, g = self.loglikelihood_and_gradient()
        
        # Now need to calculate the derivatives
        nEdge = len(self.edge_to_blen)  # number of edges
        l = nEdge / 2 + 1               # number of leaves
        k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

        other_derives, edge_derives = g[:-nEdge], g[-nEdge:]
        edge_to_derives = {self.edge_list[i] : edge_derives[i] for i in range(len(self.edge_list))}

        leaf_branch = [edge for edge in self.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
        out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
        internal_branch = [x for x in self.edge_to_blen.keys() if not x in leaf_branch]
        assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

        leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order

        Lr_derives = []  # used to store derivatives for the clock parameters L, r0, r1, ...
        Lr_derives.append(sum(edge_derives))  # dLL/dL = sum(all derives)
        Lr_derives.append(edge_to_derives[out_group_branch] * 2 / (self.x_Lr[1] - 2)
                          + sum(edge_derives))

        for i in range(2, len(self.x_Lr)):  # r(i-1)
            if self.x_Lr[i] < 1:  # when no inf could happen in the transformation
                Lr_derives.append( edge_to_derives[('N' + str(i - 2), 'N' + str(i - 1))] * (self.x_Lr[i] / (self.x_Lr[i] - 1))# 
                               + sum([edge_to_derives[internal_branch[j]] for j in range(i - 1, len(internal_branch))])  # only sum over nodes decendent from node i-1
                               + sum([edge_to_derives[leaf_branch[j]] for j in range(i - 1, len(leaf_branch))]))  # only sum over nodes decendent from node i-1
            else:  # get numerical derivative instead when inf happens
                ll = self._loglikelihood2()[0]
                self.x_clock[i + len(other_derives)] += 1e-8
                self.update_by_x_clock()
                l = self._loglikelihood2()[0]
                Lr_derives.append((l - ll) / 1e-8)
                self.x_clock[i + len(other_derives)] -= 1e-8
                self.update_by_x_clock()

        #TODO: Need to change the two sums if using general tree

        g_clock = np.concatenate( (np.array(other_derives), np.array(Lr_derives)))

        if display:
            print ('log likelihood = ', f)
            print ('Lr derivatives = ', Lr_derives)
            print ('other derivatives = ', other_derives)
            print ('Current x_clock array = ', self.x_clock)

        return f, g_clock

    def objective_wo_derivative(self, display, x):
        if self.clock:
            self.update_by_x_clock(x)
            ll = self._loglikelihood2()[0]
        else:
            self.update_by_x(x)
            ll = self._loglikelihood2()[0]

        if display:
            print ('log likelihood = ', ll)
            if self.clock:
                print ('Current x_clock array = ', self.x_clock)
            else:
                print ('Current x array = ', self.x)

        return -ll

    def objective_wo_derivative_global(self, display, x):
        if self.clock:
            self.update_by_x_clock(x, transformation = 'Exp_Neg')
            ll = self._loglikelihood2()[0]
        else:
            self.update_by_x(x, transformation = 'Exp_Neg')
            ll = self._loglikelihood2()[0]

        if display:
            print ('log likelihood = ', ll)
            if self.clock:
                print ('Current x_clock array = ', self.x_clock)
            else:
                print ('Current x array = ', self.x)

        return -ll
        
    def get_mle(self, display = True, derivative = True, em_iterations = 3, method = 'basin-hopping', niter = 2000):
        if em_iterations > 0:
            ll = self._loglikelihood2()
            # http://jsonctmctree.readthedocs.org/en/latest/examples/hky_paralog/yeast_geneconv_zero_tau/index.html#em-for-edge-lengths-only
            observation_reduction = None
            self.x_rates = np.log(optimize_em(self.get_scene(), observation_reduction, em_iterations))
            self.x = np.concatenate((self.x_process, self.x_rates))
            if self.clock:
                self.update_x_clock_by_x()
                self.update_by_x_clock()
            else:
                self.update_by_x()
                
            if display:
                print ('log-likelihood = ', ll)
                print ('updating blen length using EM')
                print ('current log-likelihood = ', self._loglikelihood2())
        else:
            if self.clock:
                self.update_by_x_clock()
            else:
                self.update_by_x()

        bnds = [(None, -0.05)] * 3
        if not self.clock:
            self.update_by_x()
            if derivative:
                f = partial(self.objective_and_gradient, display)
            else:
                f = partial(self.objective_wo_derivative, display)
            guess_x = self.x            
            bnds.extend([(None, None)] * (len(self.x_process) - 3))
            edge_bnds = [(None, None)] * len(self.x_rates)
            edge_bnds[1] = (self.minlogblen, None)
            bnds.extend(edge_bnds)
            
        else:
            self.update_by_x_clock()  # TODO: change force for blen in x_clock
            if derivative:
                f = partial(self.Clock_wrap, display)
            else:
                f = partial(self.objective_wo_derivative, display)
            guess_x = self.x_clock
            bnds.extend([(None, None)] * (len(self.x_clock) - 2 - (len(self.edge_to_blen) / 2 + 1)))
            bnds.extend([(-10, 0.0)] * (len(self.edge_to_blen) / 2))
        if method == 'BFGS':
            if derivative:
                result = scipy.optimize.minimize(f, guess_x, jac = True, method = 'L-BFGS-B', bounds = bnds)
            else:
                result = scipy.optimize.minimize(f, guess_x, jac = False, method = 'L-BFGS-B', bounds = bnds)
        elif method == 'differential_evolution':
            f = partial(self.objective_wo_derivative_global, display)
            if self.clock:
                bnds = [(np.exp(-20), 1.0 - np.exp(-20))] * len(self.x_clock)
            else:
                bnds = [(np.exp(-20), 1.0 - np.exp(-20))] * len(self.x)
                
            result = scipy.optimize.differential_evolution(f, bnds, callback = self.check_boundary_differential_evolution)
        elif method == 'basin-hopping':
            if derivative:
                result = scipy.optimize.basinhopping(f, guess_x, minimizer_kwargs = {'method':'L-BFGS-B', 'jac':True, 'bounds':bnds}, niter = niter)#, callback = self.check_boundary)
            else:
                result = scipy.optimize.basinhopping(f, guess_x, minimizer_kwargs = {'method':'L-BFGS-B', 'jac':False, 'bounds':bnds}, niter = niter)#, callback = self.check_boundary)
            
        print (result)
        self.save_x()
        return result
    def check_boundary(self, x, f, accepted):
        print("at minimum %.4f accepted %d" % (f, int(accepted)))
        return self.edge_to_blen[self.edge_list[1]] > np.exp(self.minlogblen)

    def check_boundary_differential_evolution(self, x, convergence):
        print("at lnL %.4f convergence fraction %d" % (self.ll, convergence))
        #return self.edge_to_blen[self.edge_list[1]] > np.exp(self.minlogblen)

    
    def update_x_clock_by_x(self):
        Lr = []
        x_rates = np.exp(self.x_rates)
        for i in range(len(self.edge_list) / 2 - 1):
            elist = {self.edge_list[a]:x_rates[a] for a in range(len(self.edge_list)) if self.edge_list[a][0] == 'N' + str(i)}
            elenlist = [elist[t] for t in elist]
            if i == 0:
                extra_list = [x_rates[a] for a in range(len(self.edge_list)) if self.edge_list[a][0] == 'N' + str(1) and self.edge_list[a][1][0] != 'N']
                L = (sum(elenlist) + extra_list[0]) / 2
                r0 = 2.0 - (sum(elenlist) - elist[('N0', 'N1')]) / L
                Lr.append(L)
                Lr.append(r0)
                Lr.append(extra_list[0] / (L * r0))
            else:
                Lr.append(1 - min(elenlist) / max(elenlist))
        self.x_Lr = np.array(Lr)
        self.x_clock = np.concatenate((self.x_process, np.log(self.x_Lr)))
 
    def get_geneconvTransRed(self, get_rate = False):
        row_states = []
        column_states = []
        proportions = []
        if self.Model == 'MG94':
            Qbasic = self.get_MG94Basic()
            for i, pair in enumerate(product(self.codon_nonstop, repeat = 2)):
                ca, cb = pair
                sa = self.codon_to_state[ca]
                sb = self.codon_to_state[cb]
                if ca == cb:
                    continue
                
                # (ca, cb) to (ca, ca)
                row_states.append((sa, sb))
                column_states.append((sa, sa))
                Qb = Qbasic[sb, sa]
                if isNonsynonymous(cb, ca, self.codon_table):
                    Tgeneconv = self.tau * self.omega
                else:
                    Tgeneconv = self.tau
                proportions.append(Tgeneconv / (Qb + Tgeneconv) if (Qb + Tgeneconv) >0 else 0.0)

                # (ca, cb) to (cb, cb)
                row_states.append((sa, sb))
                column_states.append((sb, sb))
                Qb = Qbasic[sa, sb]
                proportions.append(Tgeneconv / (Qb + Tgeneconv) if (Qb + Tgeneconv) >0 else 0.0)
            
        elif self.Model == 'HKY':
            Qbasic = self.get_HKYBasic()
            for i, pair in enumerate(product('ACGT', repeat = 2)):
                na, nb = pair
                sa = self.nt_to_state[na]
                sb = self.nt_to_state[nb]
                if na == nb:
                    continue

                # (na, nb) to (na, na)
                row_states.append((sa, sb))
                column_states.append((sa, sa))
                GeneconvRate = get_HKYGeneconvRate(pair, na + na, Qbasic, self.tau)
                proportions.append(self.tau / GeneconvRate if GeneconvRate > 0 else 0.0)
                

                # (na, nb) to (nb, nb)
                row_states.append((sa, sb))
                column_states.append((sb, sb))
                GeneconvRate = get_HKYGeneconvRate(pair, nb + nb, Qbasic, self.tau)
                proportions.append(self.tau / GeneconvRate if GeneconvRate > 0 else 0.0)
                
        return {'row_states' : row_states, 'column_states' : column_states, 'weights' : proportions}


    def _ExpectedNumGeneconv(self, package = 'new', display = False):
        if self.GeneconvTransRed is None:
            self.GeneconvTransRed = self.get_geneconvTransRed()

        if package == 'new':
            self.scene_ll = self.get_scene()
            requests = [{'property' : 'SDNTRAN', 'transition_reduction' : self.GeneconvTransRed}]
            j_in = {
                'scene' : self.scene_ll,
                'requests' : requests
                }        
            j_out = jsonctmctree.interface.process_json_in(j_in)

            status = j_out['status']
            ExpectedGeneconv = {self.edge_list[i] : j_out['responses'][0][i] for i in range(len(self.edge_list))}
            return ExpectedGeneconv
        else:
            print ('Need to implement this for old package')

    def _ExpectedHetDwellTime(self, package = 'new', display = False):

        if package == 'new':
            self.scene_ll = self.get_scene()
            if self.Model == 'MG94':
                heterogeneous_states = [(a, b) for (a, b) in list(product(range(len(self.codon_to_state)), repeat = 2)) if a != b]
            elif self.Model == 'HKY':
                heterogeneous_states = [(a, b) for (a, b) in list(product(range(len(self.nt_to_state)), repeat = 2)) if a != b]
            dwell_request = [dict(
                property = 'SDWDWEL',
                state_reduction = dict(
                    states = heterogeneous_states,
                    weights = [2] * len(heterogeneous_states)
                )
            )]
            
            j_in = {
                'scene' : self.scene_ll,
                'requests' : dwell_request,
                }        
            j_out = jsonctmctree.interface.process_json_in(j_in)

            status = j_out['status']
            ExpectedDwellTime = {self.edge_list[i] : j_out['responses'][0][i] for i in range(len(self.edge_list))}
            return ExpectedDwellTime
        else:
            print ('Need to implement this for old package')

    def _ExpectedHomDwellTime(self, package = 'new', display = False):

        if package == 'new':
            self.scene_ll = self.get_scene()
            if self.Model == 'MG94':
                homogeneous_states = [(a, b) for (a, b) in list(product(range(len(self.codon_to_state)), repeat = 2)) if a == b]
            elif self.Model == 'HKY':
                homogeneous_states = [(a, b) for (a, b) in list(product(range(len(self.nt_to_state)), repeat = 2)) if a == b]
            dwell_request = [dict(
                property = 'SDWDWEL',
                state_reduction = dict(
                    states = homogeneous_states,
                    weights = [2] * len(homogeneous_states)
                )
            )]
            
            j_in = {
                'scene' : self.scene_ll,
                'requests' : dwell_request,
                }        
            j_out = jsonctmctree.interface.process_json_in(j_in)

            status = j_out['status']
            ExpectedDwellTime = {self.edge_list[i] : j_out['responses'][0][i] for i in range(len(self.edge_list))}
            return ExpectedDwellTime
        else:
            print ('Need to implement this for old package')

    def _ExpectedDirectionalNumGeneconv(self, package = 'new', display = False):
        DirectionalNumGeneconvRed = self.get_directionalNumGeneconvRed()
        if package == 'new':
            self.scene_ll = self.get_scene()
            requests = [{'property' : 'SDNTRAN', 'transition_reduction' : i} for i in DirectionalNumGeneconvRed]
            assert(len(requests) == 2)  # should be exactly 2 requests
            j_in = {
                'scene' : self.scene_ll,
                'requests' : requests
                }            
            j_out = jsonctmctree.interface.process_json_in(j_in)
            status = j_out['status']
            ExpectedDirectionalNumGeneconv = {self.edge_list[i] : [j_out['responses'][j][i] for j in range(2)] for i in range(len(self.edge_list))}
            return ExpectedDirectionalNumGeneconv
        else:
            print ('Need to implement this for old package')

            
    def get_directionalNumGeneconvRed(self):
        row12_states = []
        column12_states = []
        proportions12 = []
        
        row21_states = []
        column21_states = []
        proportions21 = []
        if self.Model == 'MG94':
            Qbasic = self.get_MG94Basic()
            for i, pair in enumerate(product(self.codon_nonstop, repeat = 2)):
                ca, cb = pair
                sa = self.codon_to_state[ca]
                sb = self.codon_to_state[cb]
                if ca == cb:
                    continue
                
                # (ca, cb) to (ca, ca)
                row12_states.append((sa, sb))
                column12_states.append((sa, sa))
                Qb = Qbasic[sb, sa]
                if isNonsynonymous(cb, ca, self.codon_table):
                    Tgeneconv = self.tau * self.omega
                else:
                    Tgeneconv = self.tau
                proportions12.append(Tgeneconv / (Qb + Tgeneconv) if (Qb + Tgeneconv) >0 else 0.0)

                # (ca, cb) to (cb, cb)
                row21_states.append((sa, sb))
                column21_states.append((sb, sb))
                Qb = Qbasic[sa, sb]
                proportions21.append(Tgeneconv / (Qb + Tgeneconv) if (Qb + Tgeneconv) >0 else 0.0)
            
        elif self.Model == 'HKY':
            Qbasic = self.get_HKYBasic()
            for i, pair in enumerate(product('ACGT', repeat = 2)):
                na, nb = pair
                sa = self.nt_to_state[na]
                sb = self.nt_to_state[nb]
                if na == nb:
                    continue

                # (na, nb) to (na, na)
                row12_states.append((sa, sb))
                column12_states.append((sa, sa))
                GeneconvRate = get_HKYGeneconvRate(pair, na + na, Qbasic, self.tau)
                proportions12.append(self.tau / GeneconvRate if GeneconvRate > 0 else 0.0)
                

                # (na, nb) to (nb, nb)
                row21_states.append((sa, sb))
                column21_states.append((sb, sb))
                GeneconvRate = get_HKYGeneconvRate(pair, nb + nb, Qbasic, self.tau)
                proportions21.append(self.tau / GeneconvRate if GeneconvRate > 0 else 0.0)
                
        return [{'row_states' : row12_states, 'column_states' : column12_states, 'weights' : proportions12},
                {'row_states' : row21_states, 'column_states' : column21_states, 'weights' : proportions21}]
        
    def get_pointMutationRed(self):
        row_states = []
        col_states = []
        proportions = []
        
        if self.Model == 'MG94':
            Qbasic = self.get_MG94Basic()
            for i, pair in enumerate(product(self.codon_nonstop, repeat = 2)):
                ca, cb = pair
                sa = self.codon_to_state[ca]
                sb = self.codon_to_state[cb]
                if ca != cb:                        
                    for cc in self.codon_nonstop:
                        if cc == ca or cc == cb:
                            continue
                        sc = self.codon_to_state[cc]

                        # (ca, cb) to (ca, cc)
                        Qb = Qbasic[sb, sc]
                        if Qb != 0:
                            row_states.append((sa, sb))
                            col_states.append((sa, sc))
                            proportions.append(1.0)

                        # (ca, cb) to (cc, cb)
                        Qb = Qbasic[sa, sc]
                        if Qb != 0:
                            row_states.append((sa, sb))
                            col_states.append((sc, sb))
                            proportions.append(1.0)
                    # (ca, cb) to (ca, ca)
                    row_states.append((sa, sb))
                    col_states.append((sa, sa))
                    Qb = Qbasic[sb, sa]
                    if isNonsynonymous(cb, ca, self.codon_table):
                        Tgeneconv = self.tau * self.omega
                    else:
                        Tgeneconv = self.tau
                    proportions.append(1.0 - Tgeneconv / (Qb + Tgeneconv) if (Qb + Tgeneconv) >0 else 0.0)

                    # (ca, cb) to (cb, cb)
                    row_states.append((sa, sb))
                    col_states.append((sb, sb))
                    Qb = Qbasic[sa, sb]
                    proportions.append(1.0 - Tgeneconv / (Qb + Tgeneconv) if (Qb + Tgeneconv) >0 else 0.0)
                else:
                    for cc in self.codon_nonstop:
                        if cc == ca:
                            continue
                        sc = self.codon_to_state[cc]

                        # (ca, ca) to (ca,  cc)
                        Qb = Qbasic[sa, sc]
                        if Qb != 0:
                            row_states.append((sa, sb))
                            col_states.append((sa, sc))
                            proportions.append(1.0)
                        # (ca, ca) to (cc, ca)
                            row_states.append((sa, sb))
                            col_states.append((sc, sa))
                            proportions.append(1.0)

                        # (ca, ca) to (cc, cc)
                            row_states.append((sa, sb))
                            col_states.append((sc, sc))
                            proportions.append(1.0)
        elif self.Model == 'HKY':
            Qbasic = self.get_HKYBasic()
            for i, pair in enumerate(product('ACGT', repeat = 2)):
                na, nb = pair
                sa = self.nt_to_state[na]
                sb = self.nt_to_state[nb]
                if na!= nb:
                    for nc in 'ACGT':
                        if nc == na or nc == nb:
                            continue
                        sc = self.nt_to_state[nc]

                        # (na, nb) to (na, nc)
                        Qb = Qbasic[sb, sc]
                        if Qb != 0:
                            row_states.append((sa, sb))
                            col_states.append((sa, sc))
                            proportions.append(1.0)

                        # (na, nb) to (nc, nb)
                        Qb = Qbasic[sa, sc]
                        if Qb != 0:
                            row_states.append((sa, sb))
                            col_states.append((sc, sb))
                            proportions.append(1.0)

                    # (na, nb) to (na, na)
                    row_states.append((sa, sb))
                    col_states.append((sa, sa))
                    Qb = Qbasic[sb, sa]
                    proportions.append(Qb / (Qb + self.tau))

                    # (na, nb) to (nb, nb)
                    row_states.append((sa, sb))
                    col_states.append((sb, sb))
                    Qb = Qbasic[sa, sb]
                    proportions.append(Qb / (Qb + self.tau))
                     
        return [{'row_states' : row_states, 'column_states' : col_states, 'weights' : proportions}]

    def _ExpectedpointMutationNum(self, package = 'new', display = False):
        pointMutationRed = self.get_pointMutationRed()
        if package == 'new':
            self.scene_ll = self.get_scene()
            requests = [{'property' : 'SDNTRAN', 'transition_reduction' : i} for i in pointMutationRed]
            assert(len(requests) == 1)  # should be exactly 1 request
            j_in = {
                'scene' : self.scene_ll,
                'requests' : requests
                }            
            j_out = jsonctmctree.interface.process_json_in(j_in)
            status = j_out['status']
            ExpectedpointMutation = {self.edge_list[i] : j_out['responses'][0][i] for i in range(len(self.edge_list))}
            return ExpectedpointMutation
        else:
            print ('Need to implement this for old package')

    def _SitewiseExpectedpointMutationNum(self, package = 'new', display = False):
        pointMutationRed = self.get_pointMutationRed()
        if package == 'new':
            self.scene_ll = self.get_scene()
            requests = [{'property' : 'DDNTRAN', 'transition_reduction' : i} for i in pointMutationRed]
            assert(len(requests) == 1)  # should be exactly 1 request
            j_in = {
                'scene' : self.scene_ll,
                'requests' : requests
                }            
            j_out = jsonctmctree.interface.process_json_in(j_in)
            status = j_out['status']
            SitewiseExpectedpointMutation = np.matrix(j_out['responses'][0]).T
            return SitewiseExpectedpointMutation
        else:
            print ('Need to implement this for old package')

    def _SitewiseExpectedDirectionalNumGeneconv(self, package = 'new', display = False):
        DirectionalNumGeneconvRed = self.get_directionalNumGeneconvRed()
        if package == 'new':
            self.scene_ll = self.get_scene()
            requests = [{'property' : 'DDNTRAN', 'transition_reduction' : i} for i in DirectionalNumGeneconvRed]
            assert(len(requests) == 2)  # should be exactly 2 requests
            j_in = {
                'scene' : self.scene_ll,
                'requests' : requests
                }            
            j_out = jsonctmctree.interface.process_json_in(j_in)
            status = j_out['status']
            SitewiseExpectedDirectionalNumGeneconv = [np.matrix(j_out['responses'][0]).T, np.matrix(j_out['responses'][1]).T]
            return SitewiseExpectedDirectionalNumGeneconv
        else:
            print ('Need to implement this for old package')

    def get_SitewisePosteriorSummary(self, summary_path, file_name = None):
        SitewiseExpectedpointMutation = self._SitewiseExpectedpointMutationNum()
        SiteExpectedDirectionalNumGeneconv = self._SitewiseExpectedDirectionalNumGeneconv()
        if file_name == None:
            if not self.Force:
                prefix_summary = summary_path + self.Model + '_'
            else:
                prefix_summary = summary_path + 'Force_' + self.Model + '_'
                

            if self.clock:
                suffix_summary = '_clock_sitewise_summary.txt'
            else:
                suffix_summary = '_nonclock_sitewise_summary.txt'    

            summary_file = prefix_summary + '_'.join(self.paralog) + suffix_summary
        else:
            summary_file = file_name

        header = ' '.join(str(i + 1) for i in range(SitewiseExpectedpointMutation.shape[1]))
        label_pointMutation = ' '.join('(' + ','.join(edge) + ',mut)' for edge in self.edge_list)
        label_DirIGC = ' '.join([' '.join('(' + ','.join(edge) + ',1->2)' for edge in self.edge_list),
                        ' '.join('(' + ','.join(edge) + ',2->1)' for edge in self.edge_list)])

        np.savetxt(open(summary_file, 'w+'), np.concatenate((SitewiseExpectedpointMutation, SiteExpectedDirectionalNumGeneconv[0], SiteExpectedDirectionalNumGeneconv[1]), axis = 0),
                   delimiter = ' ', header = header, footer = ' '.join([label_pointMutation, label_DirIGC]))
            

    def get_ExpectedNumGeneconv(self):
        self.ExpectedGeneconv = self._ExpectedNumGeneconv()

    def get_ExpectedHetDwellTime(self):
        self.ExpectedDwellTime = self._ExpectedHetDwellTime()

    def numerical_Clock_derivative(self):
        ll = self._loglikelihood2()[0]
        Clock_drv = []
        for i in range(len(self.x_clock)):
            self.x_clock[i] += 1e-8
            self.update_by_x_clock()
            l = self._loglikelihood2()[0]
            Clock_drv.append((l - ll) / 1e-8)
            self.x_clock[i] -= 1e-8
            self.update_by_x_clock()
        return Clock_drv

    def get_summary(self, output_label = False):
        out = [self.nsites, self.ll]
        out.extend(self.pi)
        if self.Model == 'HKY': # HKY model doesn't have omega parameter
            out.extend([self.kappa, self.tau])
            label = ['length', 'll','pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'tau']
        elif self.Model == 'MG94':
            out.extend([self.kappa, self.omega, self.tau])
            label = ['length', 'll','pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'omega', 'tau']

        k = len(label)  # record the length of non-blen parameters

        label.extend(self.edge_list)

        out.extend([self.edge_to_blen[label[j]] for j in range(k, len(label))])

        if not self.ExpectedGeneconv:
            self.get_ExpectedNumGeneconv()

        if not self.ExpectedDwellTime:
            self.get_ExpectedHetDwellTime()

        label.extend([ (a, b, 'tau') for (a, b) in self.edge_list])
        out.extend([self.ExpectedGeneconv[i] / (self.edge_to_blen[i] * self.ExpectedDwellTime[i]) if self.ExpectedDwellTime[i] != 0 else 0 for i in self.edge_list])


        # Now add directional # of geneconv events
        ExpectedDirectionalNumGeneconv = self._ExpectedDirectionalNumGeneconv()
        label.extend([ (a, b, '1->2') for (a, b) in self.edge_list])
        out.extend([ExpectedDirectionalNumGeneconv[i][0] for i in self.edge_list])
        label.extend([ (a, b, '2->1') for (a, b) in self.edge_list])
        out.extend([ExpectedDirectionalNumGeneconv[i][1] for i in self.edge_list])


        # Now add Expected # of point mutations

        ExpectedPointMutation = self._ExpectedpointMutationNum()
        label.extend([ (a, b, 'mut') for (a, b) in self.edge_list])
        out.extend([ExpectedPointMutation[i] for i in self.edge_list])

        for i in range(k, len(label)):
            label[i] = '(' + ','.join(label[i]) + ')'

        if output_label:
            return out, label
        else:
            return out        

    def get_individual_summary(self, summary_path, file_name = None):
        if file_name == None:
            if not self.Force:
                prefix_summary = summary_path + self.Model + '_'
            else:
                prefix_summary = summary_path + 'Force_' + self.Model + '_'
                

            if self.clock:
                suffix_summary = '_clock_summary.txt'
            else:
                suffix_summary = '_nonclock_summary.txt'    

            summary_file = prefix_summary + '_'.join(self.paralog) + suffix_summary
        else:
            summary_file = file_name
        res = self.get_summary(True)
        summary = np.matrix(res[0])
        label = res[1]
            
        footer = ' '.join(label)  # row labels
        np.savetxt(open(summary_file, 'w+'), summary.T, delimiter = ' ', footer = footer)

    def gen_save_file_name(self):
        prefix_save = self.save_path + self.Model
        if self.Force:
            prefix_save = prefix_save + '_Force'

##        if self.Dir:
##            prefix_save = prefix_save + '_Dir'
##
##        if self.gBGC:
##            prefix_save = prefix_save + '_gBGC'

        if self.clock:
            suffix_save = '_clock_save.txt'
        else:
            suffix_save = '_nonclock_save.txt'

        save_file = prefix_save +'_' + '_'.join(self.paralog) + suffix_save
        return save_file

    def save_x(self):
        if self.clock:
            save = self.x_clock
        else:
            save = self.x

        if self.save_name == None:
            save_file = self.gen_save_file_name()
        else:
            save_file = self.save_name
            
        np.savetxt(open(save_file, 'w+'), save.T)

    def initialize_by_save(self, save_file):
            
        if self.clock:
            self.x_clock = np.loadtxt(open(save_file, 'r'))
            self.update_by_x_clock()
        else:
            self.x = np.loadtxt(open(save_file, 'r'))
            self.update_by_x()

def main(args):
    paralog = [args.paralog1, args.paralog2]
    #alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    alignment_file = './NewPairsAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    newicktree = '../PairsAlignemt/YeastTree.newick'
    path = './NewPackageNewRun/'
    summary_path = './NewPackageNewRun/'
    omega_guess = 0.1

    print ('Now calculate MLE for pair', paralog)
    
    if hasattr(args, 'Force'):
        Force = args.Force
        Force_hky = None
        if Force != None:
            path += 'Force_'
            Force_hky = {4:0.0}
    else:
        Force = None
        Force_hky = None

    test_hky = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = Force_hky, clock = False)
    result_hky = test_hky.get_mle(display = False, derivative = True, em_iterations = 1, method = 'BFGS')
    test_hky.get_ExpectedNumGeneconv()
    test_hky.get_ExpectedHetDwellTime()
    test_hky.get_individual_summary(summary_path = summary_path)
    test_hky.save_to_file(path = path)

    test2_hky = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = Force_hky, clock = True)
    result2_hky = test2_hky.get_mle(display = False, derivative = True, em_iterations = 1, method = 'BFGS')
    test2_hky.get_ExpectedNumGeneconv()
    test2_hky.get_ExpectedHetDwellTime()
    test2_hky.get_individual_summary(summary_path = summary_path)
    test2_hky.save_to_file(path = path)


    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = False)
    x = np.concatenate((test_hky.x_process[:-1], np.log([omega_guess]), test_hky.x_process[-1:], test_hky.x_rates))
    test.update_by_x(x)
    
    result = test.get_mle(display = True, derivative = True, em_iterations = 1, method = 'BFGS')
    test.get_ExpectedNumGeneconv()
    test.get_ExpectedHetDwellTime()
    test.get_individual_summary(summary_path = summary_path)
    test.save_to_file(path = path)

    test2 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = True)
    #x_clock = np.concatenate((test2_hky.x_process[:-1], np.log([omega_guess]), test2_hky.x_process[-1:], test2_hky.x_Lr))
    #test2.update_by_x_clock(x_clock)
    result = test2.get_mle(display = True, derivative = True, em_iterations = 0, method = 'BFGS')
    test2.get_ExpectedNumGeneconv()
    test2.get_individual_summary(summary_path = summary_path)
    test2.save_to_file(path = path)
    

if __name__ == '__main__':
##    parser = argparse.ArgumentParser()
##    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
##    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
##    parser.add_argument('--Force', type = ast.literal_eval, help = 'Parameter constraints')
##    
##    main(parser.parse_args())

##    paralog1 = 'YNL069C'
##    paralog2 = 'YIL133C'
##    paralog1 = 'YBL087C'
##    paralog2 = 'YER117W'
##    paralog1 = 'YDR502C'
##    paralog2 = 'YLR180W'
    paralog1 = 'YLR406C'
    paralog2 = 'YDL075W'
##    #paralog1 = 'YML026C'
##    #paralog2 = 'YDR450W'
####    path = './NewPackageNewRun/'
######    paralog1 = 'ECP'
######    paralog2 = 'EDN'
######
####    Force    = {5:0.0}
######
    paralog = [paralog1, paralog2]
    #alignment_file = './simulation/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '.fasta'
    alignment_file = '../../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    save_name = '../test_save.txt'
    #alignment_file = './BootStrapSamples/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_Boot1.fasta'
    #save_name = './BootStrapSave/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_Boot1_save.txt'
##    alignment_file = './TestTau/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_switched.fasta'
##    #alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
##    alignment_file = './NewPairsAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
##    #alignment_file = '../data/cleaned_input_data.fasta'
##    #alignment_file = '../data/cleanedfasta.fasta'
##    #newicktree = './YeastTree_remove_Cas.newick'
    #newicktree = '../PairsAlignemt/YeastTree.newick'
    newicktree = '../YeastTree.newick'
##    #newicktree = '../data/input_tree.newick'

##    x = np.array([-0.72980621, -0.56994663, -0.96216856,  1.73940961, -1.71054117,  0.54387332,
##                  -1.33866266, -3.47424374, -1.68831105, -1.80543811, -3.62520339, -2.69364618,
##                  -4.41935536, -2.61416486, -3.63272477, -2.78779654, -3.17653234, -3.95894103])
##
##    x_clock = np.array([-0.69711204, -0.49848498, -1.33603066,  1.861808,   -2.21507185,  5.23430255,
##                        -5.82838102, -0.07656569, -1.56484193, -0.18590612, -0.16381505, -0.41778555,
##                        -0.18178853])
####
##    out_group_blen = np.arange(0.0001, 0.01, 0.005)
##    ll_list = []

    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = None, clock = False, save_name = save_name)
    test.get_mle(False, True, 0, 'BFGS')
    print (test.tau)
    print (test.gen_save_file_name())
    print (test._loglikelihood2())
    test.get_SitewisePosteriorSummary(summary_path = '../Summary/')
#    test.get_individual_summary(summary_path = './BootStrapSummary/' + '_'.join(paralog) + '/', file_name = './BootStrapSummary/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_Boot1_summary.txt')



##    paralog1 = 'YER131W'
##    paralog2 = 'YGL189C'
##    paralog = [paralog1, paralog2]
##    alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
##
##    newicktree = '../PairsAlignemt/YeastTree.newick'
##    test3 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = None, clock = True)
##    test3.get_mle(True, True, 0, 'BFGS')
##    print test3.tau
##
##    test4 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = None, clock = True)
##    test4.get_mle(True, False, 0, 'BFGS')
##    print test4.tau

##    for blen in out_group_blen:
##        test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = {6:blen}, clock = False)
##        #test.Force = Force = {6:blen}
##        #test.update_by_x()
##        test.get_mle(False, True, 1, 'BFGS')
##        ll_list.append(test.ll)
    
##    np.savetxt(open('./testlikelihood.txt', 'w+'), np.matrix([ll_list, out_group_blen]), delimiter = ' ')
##    test.get_mle(False, False, 1, 'basin-hopping')
    
##    test3.update_by_x(x)
##    #test3.update_by_x_clock(x_clock)
##    #print test3._loglikelihood2()
##    #test3.get_mle(display = True, em_iterations = 0)
####    test3.save_to_file(path = path)
##    
##    test4 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = True)
##    test4.update_by_x_clock(x_clock)
##    test4.get_mle(display = True)
##    test4.save_to_file(path = path)
    
##    #test.get_ExpectedNumGeneconv()
##    #print test.ExpectedGeneconv
####    x = np.array([-0.64353555, -0.48264002, -0.93307917,  1.76977596, -2.44028574, 1.19617746,
####                   -3.90910406, -2.01786889,  # N0
####                   -2.92359461, -2.83555499,  # N1
####                   -4.30005794, -3.29056063,  # N2
####                   -4.59465356, -3.49891453,  # N3
####                   -6.28014701, -3.10701318,  # N4
####                   -3.87159909, -3.97438916]) # N5
####                 
####    test.update_by_x(x = x)
####    test.get_mle(display = True)
####    test.save_to_file(path = path)
####
####    # MG94 Force nonclock
####    test2 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = False, nnsites = 157)
####    test2.update_by_x(x = test.x)
####    test2.get_mle(display = True)
####    test2.save_to_file(path = path)
####
####    # MG94 clock
####    test3 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = None, clock = True, nnsites = 157)
####    test3.update_by_x_clock(x_clock = test.x_clock)
####    test3.get_mle(display = True)
####    test3.save_to_file(path = path)
####
####    # MG94 nonclock
####    test4 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = None, clock = False, nnsites = 157)
####    test4.update_by_x(x = test3.x)
####    test4.get_mle(display = True)
####    test4.save_to_file(path = path)
##
##    print test._loglikelihood()  # should be -1513.0033676643861
##    test.get_mle(False)
##    test.get_ExpectedHetDwellTime()
##    test.get_ExpectedNumGeneconv()
##
##    for i in test.edge_list:
##        if test.ExpectedDwellTime[i] != 0:
##            print i, test.ExpectedGeneconv[i]/(test.ExpectedDwellTime[i] * test.edge_to_blen[i])
##
##
####    paralog1 = 'YDR502C'
####    paralog2 = 'YLR180W'
####
####
####    Force    = {5:0.0}
####
####    paralog = [paralog1, paralog2]
####    alignment_file = '../PairsAlignemt/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
####
####    newicktree = '../PairsAlignemt/YeastTree.newick'
####    test2 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = None, clock = False)#, nnsites = 5)
####
####
####    x = np.array([-0.76374179, -0.56512517, -0.8661595 ,  1.33286957, -3.12679708, -0.55435137,
####                  -1.59080985, -1.35214307,   # N0
####                  -1.48481827, -0.5117971 ,   # N1
####                  -2.4550277 , -2.02888192,   # N2
####                  -2.47147681, -1.7627708 ,   # N3
####                  -2.62667816, -1.64821228,   # N4
####                  -2.38731535, -2.41097359])  # N5
####                  
####    test2.update_by_x(x = x)
####
####    print test2._loglikelihood()
####        
##        
##        
