# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Re-write of my previous CodonGeneconv class
# commit number: Oct 22nd, 2014 for old package
# cb1ba60ee2b57d6703cd9a3987000c2fd4dd68a5
# commit number: Dec 17th, 2014 for new package
# 33e393a973161e3a29149e82bfda23882b5826f3

  

from __future__ import print_function, absolute_import
from CodonGeneconFunc import *
import argparse
#from jsonctmctree.extras import optimize_em
import ast

#import matplotlib.pyplot as plt

class ReCodonGeneconv:
    def __init__(self, tree_newick, alignment, paralog, Model = 'MG94', nnsites = None, clock = False, Force = None, save_path = './save/', save_name = None, post_dup = 'N1'):
        self.newicktree  = tree_newick  # newick tree file loc
        self.seqloc      = alignment    # multiple sequence alignment, now need to remove gap before-hand
        self.paralog     = paralog      # parlaog list
        self.nsites      = nnsites      # number of sites in the alignment used for calculation
        self.Model       = Model
        self.IGC         = ''           # indicates if or not IGC in ancestral compare result
        self.ll          = 0.0          # Store current log-likelihood
        self.Force       = Force        # parameter constraints only works on self.x not on x_clock which should be translated into self.x first
        self.clock       = clock        # molecular clock control
        self.post_dup    = post_dup     # store first post-duplication node name
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
        self.state_to_nt    = {i:a for (i, a) in enumerate('ACGT')}
        self.codon_table    = dict(zip(codons, amino_acids))
        self.codon_nonstop  = [a for a in self.codon_table.keys() if not self.codon_table[a]=='*']
        self.codon_to_state = {a.upper() : i for (i, a) in enumerate(self.codon_nonstop)}
        self.state_to_codon = {i : a.upper() for (i, a) in enumerate(self.codon_nonstop)}
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

        # ancestral reconstruction series
        self.reconstruction_series  = None  # nodes * paralogs * 'string'

        # Initialize all parameters
        self.initialize_parameters()
        
        
    def initialize_parameters(self):
        self.get_tree()
        self.get_data()
        self.get_initial_x_process()
        save_file = self.get_save_file_name()

        if os.path.isfile(save_file):  # if the save txt file exists and not empty, then read in parameter values
            if os.stat(save_file).st_size > 0:
                self.initialize_by_save(save_file)
                print ('Successfully loaded parameter value from ' + save_file)
                
    def get_tree(self):
        self.tree, self.edge_list, self.node_to_num = read_newick(self.newicktree, self.post_dup)
        self.num_to_node = {self.node_to_num[i]:i for i in self.node_to_num}
        self.edge_to_blen = {edge:1.0 for edge in self.edge_list}


    def nts_to_codons(self):
        for name in self.name_to_seq.keys():
            assert(len(self.name_to_seq[name]) % 3 == 0)
            tmp_seq = [self.name_to_seq[name][3 * j : 3 * j + 3] for j in range(int(len(self.name_to_seq[name]) / 3) )]
            self.name_to_seq[name] = tmp_seq
       
    def get_data(self):
        seq_dict = SeqIO.to_dict(SeqIO.parse( self.seqloc, "fasta" ))
        self.name_to_seq = {name:str(seq_dict[name].seq) for name in seq_dict.keys()}
        
        if self.Model == 'MG94':
            # Convert from nucleotide sequences to codon sequences.
            self.nts_to_codons()
            obs_to_state = deepcopy(self.codon_to_state)
            obs_to_state['---'] = -1
        else:
            obs_to_state = deepcopy(self.nt_to_state)
            obs_to_state['-'] = -1

        # change the number of sites for calculation if requested
        if self.nsites is None:
            self.nsites = len(self.name_to_seq[list(self.name_to_seq.keys())[0]])
        else:
            for name in self.name_to_seq:
                self.name_to_seq[name] = self.name_to_seq[name][: self.nsites]
        print ('number of sites to be analyzed: ', self.nsites)

        # assign observable parameters
        self.observable_names = [n for n in self.name_to_seq.keys() if self.separate_species_paralog_names(n)[0] in self.node_to_num.keys()]
        suffix_to_axis = {n:i for (i, n) in enumerate(list(set(self.paralog)))}
        self.observable_nodes = [self.node_to_num[self.separate_species_paralog_names(n)[0]] for n in self.observable_names]
        self.observable_axes = [suffix_to_axis[self.separate_species_paralog_names(s)[1]] for s in self.observable_names]
        
        # Now convert alignment into state list
        iid_observations = []
        for site in range(self.nsites):
            observations = []
            for name in self.observable_names:
                observation = obs_to_state[self.name_to_seq[name][site]]
                observations.append(observation)
            iid_observations.append(observations)
        self.iid_observations = iid_observations

    def separate_species_paralog_names(self, seq_name):
        assert(seq_name in self.name_to_seq)  # check if it is a valid sequence name
        matched_paralog = [paralog for paralog in self.paralog if paralog in seq_name]
        # check if there is exactly one paralog name in the sequence name
        return [seq_name.replace(matched_paralog[0], ''), matched_paralog[0]]

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
            self.x_Lr = np.log(np.ones(int(l)) * 0.6)

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
        if not x_clock is None:
            self.x_clock = x_clock
        self.unpack_x_clock(transformation = transformation)
        self.update_by_x(transformation = transformation)
        
    def unpack_x_clock(self, transformation):
        assert(self.clock)
        nEdge = len(self.edge_to_blen)  # number of edges
        assert(nEdge % 2 == 0)
        l = int(nEdge / 2) + 1               # number of leaves
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
        if x is not None:
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
            distn = [reduce(mul, [self.pi['ACGT'.index(b)]  for b in codon], 1) for codon in self.codon_nonstop ]
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

    def _sitewise_loglikelihood(self):
        scene = self.get_scene()
        
        log_likelihood_request = {'property':'dnnlogl'}
        requests = [log_likelihood_request]
        
        j_in = {
            'scene' : self.scene_ll,
            'requests' : requests
            }
        j_out = jsonctmctree.interface.process_json_in(j_in)

        status = j_out['status']
    
        ll = j_out['responses'][0]
        self.ll = ll

        return ll

    def get_sitewise_loglikelihood_summary(self, summary_file):
        ll = self._sitewise_loglikelihood()
        with open(summary_file, 'w+') as f:
            f.write('#Site\tlnL\t\n')
            for i in range(self.nsites):
                f.write('\t'.join([str(i), str(ll[i])]) + '\n')
    
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
            print ('log likelihood = ', -f)
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
        
    def get_mle(self, display = True, derivative = True, em_iterations = 0, method = 'BFGS', niter = 2000):
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
            bnds.extend([(None, None)] * (len(self.x_process) - 4))
            bnds.extend([(None, 7.0)] * (1))  # Now add upper limit for tau
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
            assert(len(self.edge_to_blen) % 2 == 0)
            l = int(len(self.edge_to_blen) / 2)
            bnds.extend([(None, None)] * (len(self.x_clock) - 2 - ( l + 1)))
            bnds.extend([(-10, 0.0)] * l)
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

    def isSynonymous(self, first_codon, second_codon):
        return self.codon_table[first_codon] == self.codon_table[second_codon]

    def _ExpectedHetDwellTime(self, package = 'new', display = False):
        
        if package == 'new':
            self.scene_ll = self.get_scene()
            if self.Model == 'MG94':
                syn_heterogeneous_states = [(a, b) for (a, b) in list(product(range(len(self.codon_to_state)), repeat = 2)) if a != b and self.isSynonymous(self.codon_nonstop[a], self.codon_nonstop[b])]
                nonsyn_heterogeneous_states = [(a, b) for (a, b) in list(product(range(len(self.codon_to_state)), repeat = 2)) if a != b and not self.isSynonymous(self.codon_nonstop[a], self.codon_nonstop[b])]
                dwell_request = [dict(
                    property='SDWDWEL',
                    state_reduction=dict(
                        states=syn_heterogeneous_states,
                        weights=[2] * len(syn_heterogeneous_states)
                    )),
                    dict(
                        property='SDWDWEL',
                        state_reduction=dict(
                            states=nonsyn_heterogeneous_states,
                            weights=[2] * len(nonsyn_heterogeneous_states)
                        ))
                ]

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

            ExpectedDwellTime = [{self.edge_list[i] : j_out['responses'][j][i] for i in range(len(self.edge_list))} for j in range(len(j_out))]
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
                else:
                    for nc in 'ACGT':
                        if nc == na:
                            continue
                        sc = self.nt_to_state[nc]

                        Qb = Qbasic[sa, sc]
                        if Qb != 0.0:
                            row_states.append((sa, sb))
                            col_states.append((sa, sc))
                            proportions.append(1.0)

                            row_states.append((sa, sb))
                            col_states.append((sc, sa))
                            proportions.append(1.0)                           

                     
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

        self._get_dwell_time_summary(out, label)

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

    def _get_dwell_time_summary(self, out, label):
        if not self.ExpectedGeneconv:
            self.get_ExpectedNumGeneconv()

        if not self.ExpectedDwellTime:
            self.get_ExpectedHetDwellTime()

        if self.Model == 'HKY':
            label.extend([ (a, b, 'tau') for (a, b) in self.edge_list])
            out.extend([self.ExpectedGeneconv[i] / (self.edge_to_blen[i] * self.ExpectedDwellTime[0][i]) if self.ExpectedDwellTime[0][i] != 0 else 0 for i in self.edge_list])
        elif self.Model == 'MG94':
            label.extend([(a, b, 'syn_dwell') for (a, b) in self.edge_list])
            out.extend([self.edge_to_blen[i] * self.ExpectedDwellTime[0][i] if
                        self.ExpectedDwellTime[0][i] != 0 else 0 for i in self.edge_list])
            label.extend([(a, b, 'nonsyn_dwell') for (a, b) in self.edge_list])
            out.extend([self.edge_to_blen[i] * self.ExpectedDwellTime[1][i] if
                        self.ExpectedDwellTime[1][i] != 0 else 0 for i in self.edge_list])
            label.extend([(a, b, 'tau') for (a, b) in self.edge_list])
            out.extend([self.ExpectedGeneconv[i] / (self.edge_to_blen[i] * (self.ExpectedDwellTime[0][i]+self.omega*self.ExpectedDwellTime[1][i])) \
                            if (self.ExpectedDwellTime[0][i]+self.omega*self.ExpectedDwellTime[1][i]) != 0 else 0 for i in self.edge_list])
            label.extend([(a, b, 'test_tau') for (a, b) in self.edge_list])
            out.extend([self.ExpectedGeneconv[i] / (self.edge_to_blen[i] * (self.ExpectedDwellTime[0][i]+self.ExpectedDwellTime[1][i])) \
                            if (self.ExpectedDwellTime[0][i]+self.ExpectedDwellTime[1][i]) != 0 else 0 for i in self.edge_list])

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

    def get_save_file_name(self):
        if self.save_name is None:
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
        else:
            save_file = self.save_name
        return save_file

    def save_x(self):
        if self.clock:
            save = self.x_clock
        else:
            save = self.x

        save_file = self.get_save_file_name()
            
        np.savetxt(save_file, save.T)

    def initialize_by_save(self, save_file):
            
        if self.clock:
            self.x_clock = np.loadtxt(open(save_file, 'r'))
            self.update_by_x_clock()
        else:
            self.x = np.loadtxt(open(save_file, 'r'))
            self.update_by_x()     
    
if __name__ == '__main__':
    paralog = ['YLR406C', 'YDL075W']
    Force = None
    alignment_file = '../test/YLR406C_YDL075W_test_input.fasta'
    newicktree = '../test/YeastTree.newick'
    ##    test.get_individual_summary(summary_path = '../test/Summary/')
    ##    test.get_SitewisePosteriorSummary(summary_path = '../test/Summary/')
    # Force MG94:{5:0.0} HKY:{4:0.0}

    #MG94+tau
    MG94_tau = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = False, save_path = '../test/save/')
    # MG94_tau.get_mle(True, True, 0, 'BFGS')
    MG94_tau._loglikelihood2()
    MG94_tau.get_summary(True)
    # MG94_tau.site_reconstruction()
    # MG94_tau_series = MG94_tau.reconstruction_series
##    
##    #MG94
##    MG94 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = {5:0.0}, clock = None, save_path = '../test/save/')
##    MG94.get_mle(True, True, 0, 'BFGS')
##    MG94.site_reconstruction()
##    MG94_series = MG94.reconstruction_series
##    result = MG94_tau.find_differences_between(MG94_tau_series, MG94_series)
##    print(result)
##    
##    #HKY+tau
##    HKY_tau = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = Force, clock = None, save_path = '../test/save/')
##    MG94_tau.get_mle(True, True, 0, 'BFGS')
##    HKY_tau.site_reconstruction()
##    HKY_tau_series = HKY_tau.reconstruction_series
##    
##    #MG94
##    HKY = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = {4:0.0}, clock = None, save_path = '../test/save/')
##    MG94.get_mle(True, True, 0, 'BFGS')
##    HKY.site_reconstruction()
##    HKY_series = HKY.reconstruction_series
##    result = HKY_tau.find_differences_between(HKY_tau_series, HKY_series)
##    print(result)
    
######################################################################################
######################################################################################
######################################################################################
    
    # paralog = ['EDN', 'ECP']
    # Force = None
    # ##    alignment_file = '../test/EDN_ECP_Cleaned.fasta'
    # ##    newicktree = '../test/input_tree.newick'
    # alignment_file = '../test/EDN_ECP_Outgroup_test.fasta'
    # newicktree = '../test/Outgroup_test_tree.newick'
    # Force = None
    # model = 'MG94'
    # save_name = '../test/save/Ind_' + model + '_EDN_ECP_nonclock_Outgroup_test_save.txt'
    #
    # test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = None, save_path = '../test/save/', save_name = save_name)
    # test.get_mle()
##    test.get_mle(True, True, 0, 'BFGS')
##    test.get_individual_summary(summary_path = '../test/Summary/')
##    test.get_SitewisePosteriorSummary(summary_path = '../test/Summary/')
    # test._loglikelihood2()
    #scene = test.get_scene()
    #test.update_by_x(np.concatenate((np.log([0.1, 0.9, 0.3, 11.0, 3.4]), test.x_rates)))

    
