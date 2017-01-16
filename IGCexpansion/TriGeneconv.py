# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Re-write of my previous CodonGeneconv class
# commit number: Oct 22nd, 2014 for old package
# cb1ba60ee2b57d6703cd9a3987000c2fd4dd68a5
# commit number: Dec 17th, 2014 for new package
# 33e393a973161e3a29149e82bfda23882b5826f3

from TriGeneconvFunc import *
import argparse

class TriGeneconv:
    def __init__(self, tree_newick, alignment, 
                 save_path,  # new features
                 oldest_paralog,  # Ancient paralog
                 paralog = ['ADH1A', 'ADH1B', 'ADH1C'], Model = 'HKY', nnsites = None,
                 Force = None, clock = None, Dir = None, gBGC = None, Dis = 'Free'):
        self.newicktree     = tree_newick     # newick tree file loc
        self.seqloc         = alignment       # multiple sequence alignment, gaps are not allowed in this version
        self.paralog        = paralog         # paralogs, order matters
        self.oldest_paralog = oldest_paralog  # Ancient paralog which give birth to the first duplication event
        self.nsites         = nnsites         # number of sites in the alignment used for calculation
        self.Model          = Model           # Basic point mutation model, now only HKY85 is considered
        self.ll             = 0.0             # Store current log-likelihood
        self.Force          = Force           # Force tau zero switch
        self.clock          = clock           # strict molecular clock switch
        self.Dir            = Dir             # directional model switch
        self.gBGC           = gBGC            # GC-biased model switch
        self.save_path       = save_path        # Txt file used to save current status, legendary save & load method


        ##Not sure if I want to keep these controls, comment out but not delete
##        self.logzero     = -15.0        # used to avoid log(0), replace log(0) with -15
##        self.infinity    = 1e6          # used to avoid -inf in gradiance calculation of the clock case
##        self.minlogblen  = -9.0         # log value, used for bond constraint on branch length estimates in get_mle() function

        # Tree topology related variable
        # The tree for this project is fixed, but keep the read-in feature anyway
        self.tree         = None        # store the tree dictionary used for json likelihood package parsing
        self.edge_to_blen = None        # dictionary store the unpacked tree branch length information {(node_from, node_to):blen}
        self.edge_list    = None        # kept all edges in the same order with x_rates
        self.node_to_num  = None        # dictionary used for translating tree info from self.edge_to_blen to self.tree
        self.num_to_node  = None        # dictionary used for translating tree info from self.tree to self.edge_to_blen

        # Constants for Sequence operations
        bases = 'tcag'.upper()
        self.nt_to_state    = {a:i for (i, a) in enumerate('ACGT')}
        self.triple_to_state  = {triple:i for i, triple in enumerate(product('ACGT', repeat = 3))}

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
        self.tau            = None      # list of 6 real values, in order of
                                        # 1 -> 2, 1 -> 3, 2 -> 3, 2 -> 1, 3 -> 1, 3 -> 2, oldest -> dup, dup -> oldest
        self.num_free_tau   = None      # used to store # of tau's in x_process
        self.num_other      = None      # used to store # of other parameters in x_process
        self.Dis            = Dis       # Distance model: 'Free' or 'exponential'
        self.delta          = None      # distance parameter for exponential decay
        self.gamma          = None      # gBGC parameter


        self.processes      = None      # list of basic and geneconv rate matrices. Each matrix is a dictionary used for json parsing
                                        # basic para + tau + delta + gamma

        self.scene_ll       = None      # used for lnL calculation
        self.auto_save      = 0         # used for auto-save

        # Prior distribution on the root
        self.prior_feasible_states  = None
        self.prior_distribution     = None


##        # Expected_Geneconv Events
##        # Not sure how to deal with them yet
##        self.GeneconvTransRed  = None    # dictionary of Geneconv transition matrix used for json parsing
##        self.ExpectedGeneconv  = None    # dictionary storing expected number of geneconv events on each branch
##        self.ExpectedDwellTime = None    # dictionary storing expected total dwell time of heterogeneous states of each branch same order as self.edge_list

        # Initialize all parameters
        self.initialize_parameters()


    def initialize_parameters(self):
        self.get_tree()  # Assume the tree file stays unchanged which should be.
        self.get_data()  # Multiple sequence alignment should also stay the same.
        save_file = self.gen_save_file_name()
        if os.path.isfile(save_file):  # if the save txt file exists and not empty, then read in parameter values
            if os.stat(save_file).st_size > 0:                
                self.initialize_by_save()
                print 'Successfully loaded parameter value from ' + save_file
        else:
            self.get_initial_x_process()   # otherwise, initialize all parameter to default initial values


    def get_tree(self):
        tree = Phylo.read( self.newicktree, "newick")
        #set node number for nonterminal nodes and specify root node
        numInternalNode = 2  # Add two more nodes at root for 2 duplication events
        for clade in tree.get_nonterminals():
            clade.name = 'N' + str(numInternalNode)
            numInternalNode += 1
        tree_phy = tree.as_phyloxml(rooted = 'True')
        tree_nx = Phylo.to_networkx(tree_phy)

        triples = [(u.name, v.name, d['weight']) for (u, v, d) in tree_nx.edges(data = True)] # data = True to have the blen as 'weight'
        triples.extend([('N0', 'N1', d['weight']), ('N1', 'N2', d['weight'])])
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
        l = nEdge / 2                   # number of leaves
        k = l + 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

        leaf_branch = [edge for edge in self.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
        # Seems out_group_branch is never used
        # out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
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
        tree_process = [1 if e[0] == 'N0' and e[1] == 'N1' else 2 for e in edge_list]
        self.edge_list = edge_list

        self.tree = dict(
            row = tree_row,
            col = tree_col,
            process = tree_process,
            rate = np.ones(len(tree_row))
            )    

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
        print 'number of sites to be analyzed: ', self.nsites

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

    def assign_num_free_tau(self):
        num_tau = 1   # default one value tau model
        num_other = 0
        if self.Dis == 'Free':
            num_tau = 4
        elif self.Dis == 'exponential':
            assert(self.delta)  # self.delta should be a real value for the exponential decay
            num_other += 1  # add self.delta
        elif self.Dis == 'None':
            print 'Distance parameter not included in the model'
        else:
            print 'Check Distance Model pleaze'
            assert(0==1)

        if self.gBGC:
            assert(self.gamma == None)
            num_other += 1  # add self.gamma counting for gBGC

        if self.Dir:
            num_tau = num_tau * 2

        self.num_free_tau = num_tau
        self.num_other    = num_other
    
    def check_x_process(self):
        # self.tau list directions
        # 1 -> 2, 1 -> 3, 2 -> 3, 2 -> 1, 3 -> 1, 3 -> 2
        if self.Model == 'HKY':
            num_basic_para = 4  # %AG, %A, %C, kappa
        else:
            num_basic_para = 0
        
        total_x_process_len = num_basic_para + self.num_free_tau + self.num_other
        assert(len(self.x_process) == total_x_process_len)

        # Now check force tau cases
        if self.Force:
            assert(sum(self.tau) == 0.0)

    # add gamma parameter which is equivalent to B in Lartillot's 2013 MBE paper (Phylogenetic Patterns of GC-Biased Gene Conversion...)
    # gamma parameter isn't log-transformed in self.x but used itself
    def get_initial_x_process(self):
        count = np.array([0, 0, 0, 0], dtype = float) # count for A, C, G, T in all seq
        for name in self.name_to_seq:
            for i in range(4):
                count[i] += ''.join(self.name_to_seq[name]).count('ACGT'[i])
        count = count / count.sum()

        self.assign_num_free_tau()
        if self.Model == 'HKY':
            # deal with basic model parameters first
            self.omega = None
            self.x_process = np.log(np.array([count[0] + count[2], count[0] / (count[0] + count[2]), count[1] / (count[1] + count[3]),self.kappa]))

            # start to add in IGC specific parameters
            self.x_process = np.concatenate((self.x_process, np.log([0.1] * self.num_free_tau)))  # tau is list of num_tau elements
            if self.Dis == 'exponential':  # add in tau + delta
                self.delta = 0.1
                self.x_process = np.concatenate((self.x_process, [np.log(self.delta)]))
        
            if self.gBGC:
                self.gamma = 1.0
                self.x_process = np.concatenate((self.x_process, [self.gamma]))  # add in gamma for gBGC
        if self.Force:
            self.assign_force()
            self.tau = [0.0] * 8
        self.check_x_process()
        
        self.x_rates = np.log(np.array([ 0.01 * self.edge_to_blen[edge] for edge in self.edge_to_blen.keys()]))
        self.x = np.concatenate((self.x_process, self.x_rates))

##        if self.clock:   # set-up x_clock if it's a clock model
##            l = len(self.edge_to_blen) / 2 + 1               # number of leaves
##            self.x_Lr = np.log(np.ones((l)) * 0.9)
##
##            if transformation == 'log':
##                self.x_clock = np.concatenate((self.x_process, self.x_Lr))
##            elif transformation == 'None':
##                self.x_Lr = np.exp(self.x_Lr)
##            elif transformation == 'Exp_Neg':
##                self.x_Lr = np.exp(-np.exp(self.x_Lr))
##                
##            self.x_clock = np.concatenate((self.x_process, self.x_Lr))
##            self.unpack_x_clock(transformation = transformation)
##
##        self.update_by_x(transformation = transformation)


    def update_by_x(self, x = None):
        k = len(self.edge_to_blen)
        if not x == None:
            self.x = x
        self.x_process, self.x_rates = self.x[:-k], self.x[-k:]
        Force_process = None
        Force_rates = None
        if self.Force:
            self.assign_force()
            Force_process = {i:self.Force[i] for i in self.Force.keys() if i < len(self.x) - k}
            Force_rates = {(i-len(self.x_process)):self.Force[i] for i in self.Force.keys() if not i < len(self.x) - k}

        self.unpack_x_process(Force_process = Force_process)
        self.unpack_x_rates(Force_rates = Force_rates)

    def assign_force(self):
        # assign Force dictionary which acts on Tau based on the model setup
        # x_process = [%AG, %A, %C, kappa] + [tau] + [delta, gamma]
        if self.Model == 'HKY':
            self.Force = {(4 + i):0.0 for i in range(self.num_free_tau)}

    def unpack_x_process(self, Force_process = None):
        if self.gBGC:
            x_process = np.concatenate((np.exp(self.x_process[:-1]), [self.x_process[-1]]))
        else:
            x_process = np.exp(self.x_process)

        if Force_process != None:
            for i in Force_process.keys():
                x_process[i] = Force_process[i]

        if self.Model == 'HKY':
            # x_process = [%AG, %A, %C, kappa] + [tau] + [delta, gamma]
            pi_a = x_process[0] * x_process[1]
            pi_c = (1 - x_process[0]) * x_process[2]
            pi_g = x_process[0] * (1 - x_process[1])
            pi_t = (1 - x_process[0]) * (1 - x_process[2])
            self.pi = [pi_a, pi_c, pi_g, pi_t]
            self.kappa = x_process[3]
            self.tau = self.assign_tau(x_process[4:(4 + self.num_free_tau)])

            if self.Dis == 'exponential':
                self.delta = x_process[4 + self.num_free_tau]

            if self.gBGC:  # gamma is always stored at last position
                self.gamma = x_process[-1]

            self.check_x_process()
        # Now update the prior distribution
        self.get_prior()

        # Now update processes (Rate matrices)
        self.get_processes()        

    def assign_tau(self, x_tau): # x_tau is real value, no need to exp()
        if self.num_free_tau == 1:
            return [x_tau] * 8
        elif self.num_free_tau == 4:
            # 1 -> 2, 1 -> 3, 2 -> 3, 2 -> 1, 3 -> 1, 3 -> 2, oldest -> dup, dup -> oldest
            #tau_list = [val for val in x_tau[:-1] for _ in (0,1)]
            tau_list = [x_tau[0], x_tau[2], x_tau[1], x_tau[0], x_tau[2], x_tau[1]]
            tau_list.extend([x_tau[-1]] * 2)
            return tau_list
        elif self.num_free_tau == 8:
            return x_tau
        else:
            print "This case hasn't been implemented, please change assign_tau() function "
            assert(0==1)

    def unpack_x_rates(self, Force_rates = None):  # TODO: Change it to fit general tree structure rather than cherry tree
        x_rates = np.exp(self.x_rates)

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

    def get_prior(self):
        if self.Model == 'HKY':
            self.prior_feasible_states = [(self.nt_to_state[nt], self.nt_to_state[nt], self.nt_to_state[nt]) for nt in 'ACGT']
            distn = [ self.pi['ACGT'.index(nt)] for nt in 'ACGT' ]
            distn = np.array(distn) / sum(distn)
        self.prior_distribution = distn

    def get_scene(self):
        if self.Model == 'MG94':
            state_space_shape = [61, 61, 61]
        elif self.Model == 'HKY':
            state_space_shape = [4, 4, 4]
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

    def get_processes(self):
        #print self.tau
        oldest_paralog_num = self.paralog.index(self.oldest_paralog)
        if self.Model == 'HKY':
            self.processes = get_HKYGeneconv(self.pi, self.kappa, self.prior_distribution, self.nt_to_state, self.tau, oldest_paralog_num)

    def _loglikelihood(self, store = True, edge_derivative = False):
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

    def loglikelihood_and_gradient(self, display = False):
        '''
        Modified from Alex's objective_and_gradient function in ctmcaas/adv-log-likelihoods/mle_geneconv_common.py
        '''
        self.update_by_x()
        delta = 1e-8
        x = deepcopy(self.x)  # store the current x array

        fn = self._loglikelihood

        ll, edge_derivs = fn(edge_derivative = True)
        
        m = len(self.x) - len(self.edge_to_blen)

        # use finite differences to estimate derivatives with respect to these parameters
        other_derivs = []
        
        for i in range(m):
            if self.Force:
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
            print 'log likelihood = ', ll
            print 'Edge derivatives = ', edge_derivs
            print 'other derivatives:', other_derivs
            print 'Current x array = ', self.x

        self.ll = ll
        f = -ll
        g = -np.concatenate((other_derivs, edge_derivs))
        return f, g

    def objective_with_gradient(self, display, x):
        self.update_by_x(x)
        f, g = self.loglikelihood_and_gradient(display = display)
        self.auto_save += 1
        if self.Model == 'HKY' and self.auto_save == 5:
            self.save_x()
            self.auto_save = 0
        return f, g

    def get_mle(self, display = True, derivative = True, method = 'BFGS', niter = 2000):
        if self.clock:
            self.update_by_x_clock()
        else:
            self.update_by_x()

        bnds = [(None, -0.05)] * 3
        if not self.clock:
            self.update_by_x()
            if derivative:
                f = partial(self.objective_with_gradient, display)
##            else:
##                f = partial(self.objective_wo_derivative, display)
            guess_x = self.x            
            bnds.extend([(None, None)] * (len(self.x_process) - 3))
            edge_bnds = [(None, None)] * len(self.x_rates)
            #edge_bnds[1] = (-9.0, None)  # minlogblen
            bnds.extend(edge_bnds)
            
##        else:
##            self.update_by_x_clock()  # TODO: change force for blen in x_clock
##            if derivative:
##                f = partial(self.Clock_wrap, display)
##            else:
##                f = partial(self.objective_wo_derivative, display)
##            guess_x = self.x_clock
##            bnds.extend([(None, None)] * (len(self.x_clock) - 2 - (len(self.edge_to_blen) / 2 + 1)))
##            bnds.extend([(-10, 0.0)] * (len(self.edge_to_blen) / 2))
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

    def gen_save_file_name(self):
        prefix_save = self.save_path + self.Model
        if self.Force:
            prefix_save = prefix_save + '_Force'

        if self.Dir:
            prefix_save = prefix_save + '_Dir'

        if self.gBGC:
            prefix_save = prefix_save + '_gBGC'

        if self.clock:
            suffix_save = '_clock_save.txt'
        else:
            suffix_save = '_nonclock_save.txt'

        save_file = prefix_save + '_oldest_' + self.oldest_paralog + '_Dis_' + self.Dis + suffix_save
        return save_file

    def save_x(self):
        if self.clock:
            save = self.x_clock
        else:
            save = self.x
        np.savetxt(open(self.gen_save_file_name(), 'w+'), save.T)

    def initialize_by_save(self):
        save_file = self.gen_save_file_name()
        self.assign_num_free_tau()
        if self.clock:
            self.x_clock = np.loadtxt(open(save_file, 'r'))
            self.update_by_x_clock()
        else:
            self.x = np.loadtxt(open(save_file, 'r'))
            self.update_by_x()

    def get_summary(self, output_label = False):
        out = [self.nsites, self.ll]
        out.extend(self.pi)

        if self.Model == 'HKY': # HKY model doesn't have omega parameter
            out.extend([self.kappa])
            # 1 -> 2, 1 -> 3, 2 -> 3, 2 -> 1, 3 -> 1, 3 -> 2, oldest -> dup, dup -> oldest
            label = ['length', 'll','pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'tau12', 'tau13', 'tau23', 'tau21', 'tau31', 'tau32', 'tauOD', 'tauDO']
            
        out.extend(self.tau)
        if self.gBGC:
            out.extend([self.gamma])
            label.append('gBGC')
            
        k = len(label)  # record the length of non-blen parameters

        label.extend(self.edge_list)

        out.extend([self.edge_to_blen[label[j]] for j in range(k, len(label))])

        for i in range(k, len(label)):
            label[i] = '(' + ','.join(label[i]) + ')'

        if output_label:
            return out, label
        else:
            return out

    def get_individual_summary(self):

        summary_file = self.gen_save_file_name().replace('save', 'summary')   

        res = self.get_summary(True)
        summary = np.matrix(res[0])
        label = res[1]
            
        footer = ' '.join(label)  # row labels
        np.savetxt(open(summary_file, 'w+'), summary.T, delimiter = ' ', footer = footer)
        
if __name__ == '__main__':

    alignment_file = '../../TriGeneconv/data/ADH_intron_input.fasta'
    newicktree = '../../TriGeneconv/data/ADH1Class_tree.newick'
    save_path = '../../TriGeneconv/save/'
    paralog = ['ADH1A', 'ADH1B', 'ADH1C']
    #paralog = [paralog[0]]
    oldest_paralog = 'ADH1C'
    Dis = 'None'
    Dir = False
    gBGC = False
    Force = False
    print(oldest_paralog)

    #for oldest_paralog in paralog:
    test = TriGeneconv( newicktree, alignment_file, save_path, oldest_paralog, Force = False, Dis = Dis, Dir = Dir, gBGC = gBGC)
    test.update_by_x()
    print test._loglikelihood()
    test.get_mle(True, True)
    test.get_individual_summary()

    
    self = test







