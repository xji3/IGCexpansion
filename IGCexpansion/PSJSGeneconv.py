# A separate file for JSLikelihood
# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Xiang Ji
# xji3@ncsu.edu
from __future__ import print_function
from itertools import product
import jsonctmctree.ll, jsonctmctree.interface
from Data import Data
from Tree import Tree
from PSJSModel import PSJSModel
from Func import *
from copy import deepcopy
from functools import partial
import scipy
import scipy.optimize
import os
from Common import *


from TriGeneconv import *

class PSJSGeneconv:
    auto_save_step = 2
    def __init__(self, alignment_file, gene_to_orlg_file, # Data input
                 seq_index_file, cdna, allow_same_codon,  # Data input                 
                 tree_newick, DupLosList,                 # Tree input
                 x_js, pm_model, IGC_pm, rate_variation,  # PSJSModel input
                 node_to_pos, terminal_node_list,         # Configuration input
                 save_file, log_file,                     # Auto save files
#                root_by_dup = False,                     # JSModel input
                 force = None,                            # PSJS Parameter value constraint
                 nsites = None, space_list = None
                 ):  

        
        self.tree = Tree(tree_newick, DupLosList, terminal_node_list, node_to_pos)
        self.data = Data(alignment_file, gene_to_orlg_file, seq_index_file = seq_index_file, two_sites = True, space_list = space_list,
                         cdna = cdna, allow_same_codon = allow_same_codon)        
        self.psjsmodel = PSJSModel(x_js, pm_model, self.tree.n_orlg, IGC_pm, rate_variation, force)
        self.node_to_pos = node_to_pos
        self.terminal_node_list = terminal_node_list
        self.root_by_dup = self.tree.is_duplication_node(self.tree.phylo_tree.root.name)
        self.save_file = save_file
        self.log_file  = log_file
        self.force = force
        self.x     = None
        if os.path.isfile(self.save_file):
            self.initialize_by_save()
            print ('Loaded paramters from ' + self.save_file)
        else:            
            if self.root_by_dup:
                self.x = np.concatenate((x_js, np.log([0.01] * len(self.tree.edge_list))))
            else:
                self.x = np.concatenate((x_js, np.log([0.01] * (len(self.tree.edge_list) - 1))))

        assert(self.psjsmodel.n_orlg == self.tree.n_orlg)  # make sure n_orlg got updated
        self.ll   = None              # used to store log likelihood
        self.ExpectedGeneconv = None  # Used to store expectedNumGeneconv
        self.auto_save = 0            # used to control auto save frequency

        self.nsites = nsites
        self.iid_observations = None  # Store iid_observations which only needs one time calculation
        self.observable_nodes = None
        self.observable_axes  = None

        assert(self.self_check())


    def self_check(self):
        check_status = True
        if self.psjsmodel.rate_variation:
            check_status = check_status and self.data.cdna

        # See if the log file is there
        # It's not really a log but store values in each iteration
        if not os.path.isfile(self.log_file):
            with open(self.log_file, 'w+') as f:
                f.write('\t'.join(['#lnL ', 'x[] ', ' df ']) + '\n')
            print('Created log file: ' + self.log_file)
        else:
            print('Found log file at ' + self.log_file)
            
        return check_status

            
    def unpack_x(self, x):
        assert(len(x) == len(self.x)) # length must match
        self.x = x
        if self.root_by_dup:
            x_rate = x[-len(self.tree.edge_list):]
            x_js   = x[:(-len(self.tree.edge_list))]
        else:
            x_rate = x[-(len(self.tree.edge_list) - 1):]
            x_js   = x[:-(len(self.tree.edge_list) - 1)]
        self.unpack_x_rates(x_rate)
        self.psjsmodel.update_by_x_js(x_js)

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

    def get_prior(self):
        prior_feasible_states = []
        distn = []
        for nt_a in range(4): # Only consider nucleotide model
            for nt_b in range(4):
                js_state = (nt_a, nt_b, nt_a, nt_b)
                prior_feasible_states.append(translate_four_nt_to_two_state(js_state))
                distn.append(self.psjsmodel.PMModel.get_stationary_distn(nt_a) * self.psjsmodel.PMModel.get_stationary_distn(nt_b))

        distn = np.array(distn) / sum(distn)
        return prior_feasible_states, distn
    
    def get_scene(self, n, codon_site_pair = None):
        if self.psjsmodel.rate_variation:
            assert(not codon_site_pair == None)
        state_space_shape = self.psjsmodel.state_space_shape
        conf_list = count_process(self.tree.node_to_conf)
        if codon_site_pair == None:
            process_definitions = [self.psjsmodel.get_PM_process_definition(), self.psjsmodel.get_IGC_process_definition(n)]
        else:
            process_definitions = [self.psjsmodel.get_PM_process_definition(codon_site_pair), self.psjsmodel.get_IGC_process_definition(n, codon_site_pair)]
        self.tree.get_tree_process(conf_list)

        prior_feasible_states, distn = self.get_prior()
        
        if self.iid_observations == None:
            self.cal_iid_observations()

        if self.psjsmodel.rate_variation:
            iid_observations = self.iid_observations[codon_site_pair][n]
        else:
            assert(codon_site_pair == (1, 1) or codon_site_pair == None)
            iid_observations = self.iid_observations[n]
            
        scene = dict(
            node_count = len(self.tree.node_to_num),
            process_count = len(process_definitions),
            state_space_shape = state_space_shape,
            tree = self.tree.tree_json,
            root_prior = {'states':prior_feasible_states,
                          'probabilities':distn},
            process_definitions = process_definitions,
            observed_data = {
                'nodes':self.observable_nodes,
                'variables':self.observable_axes,
                'iid_observations':iid_observations
                }
            )
        return scene

    def cal_iid_observations(self):
        if self.iid_observations == None:
            if self.data.cdna:
                self.observable_nodes, self.observable_axes, self.iid_observations = get_all_PS_iid_observations(self.data, self.tree, self.psjsmodel.PMModel.data_type)
                if not self.psjsmodel.rate_variation:
                    new_iid_observations = {n:[] for n in self.data.space_list}
                    for n in self.data.space_list:
                        new_iid_observations[n] = [i for j in product(range(1, 4), repeat = 2) if n in self.iid_observations[j] for i in self.iid_observations[j][n] ]
                    self.iid_observations = new_iid_observations
            else:
                if self.nsites is None:
                    self.observable_nodes, self.observable_axes, self.iid_observations = get_all_PS_iid_observations(self.data, self.tree, self.psjsmodel.PMModel.data_type)
                else:
                    self.observable_nodes, self.observable_axes, self.iid_observations = get_all_PS_iid_observations(self.data, self.tree, self.psjsmodel.PMModel.data_type, self.nsites)



    def _loglikelihood(self, n, codon_site_pair, edge_derivative = False):
        if self.psjsmodel.rate_variation:
            assert(not codon_site_pair == None)
            scene = self.get_scene(n, codon_site_pair)
        else:
            scene = self.get_scene(n)
            
        log_likelihood_request = {'property':'snnlogl'}
        derivatives_request = {'property':'sdnderi'}
        if edge_derivative:
            requests = [log_likelihood_request, derivatives_request]
        else:
            requests = [log_likelihood_request]
        j_in = {
            'scene' : scene,
            'requests' : requests
            }
        j_out = jsonctmctree.interface.process_json_in(j_in)

        status = j_out['status']
    
        ll = j_out['responses'][0]
        if edge_derivative:
            edge_derivs = j_out['responses'][1]
        else:
            edge_derivs = []

        return ll, edge_derivs

    def _loglikelihood_for_one_n(self, n, edge_derivative = False):
        assert(n in self.data.space_list)
        if self.psjsmodel.rate_variation:
            ll = 0.0
            edge_derivs = 0.0
            for codon_site_pair in self.data.space_codon_site_pair[n]:
                new_ll, new_edge_derivs = self._loglikelihood(n, codon_site_pair, edge_derivative)
                ll += new_ll
                edge_derivs += np.array(new_edge_derivs)
        else:
            ll, edge_derivs = self._loglikelihood(n, codon_site_pair = None, edge_derivative = edge_derivative)

        return ll, edge_derivs

    def _sitewise_loglikelihood(self, n, codon_site_pair):
        if self.psjsmodel.rate_variation:
            assert(not codon_site_pair == None)
            scene = self.get_scene(n, codon_site_pair)
        else:
            scene = self.get_scene(n)
            
        log_likelihood_request = {'property':'dnnlogl'}
        requests = [log_likelihood_request]

        j_in = {
            'scene' : scene,
            'requests' : requests
            }
        j_out = jsonctmctree.interface.process_json_in(j_in)

        status = j_out['status']
        ll = j_out['responses'][0]

        if self.psjsmodel.rate_variation:
            return ll, self.data.space_idx_pairs[codon_site_pair][n]
        else:
            idx_pairs = [i for j in product(range(1, 4), repeat = 2) if n in self.data.two_sites_name_to_seq[j] for i in self.data.space_idx_pairs[j][n] ]
            return ll, idx_pairs

    def _sitewise_loglikelihood_for_all_n(self):
        pair_to_lnL = dict()
        for n in self.data.space_list:
            if self.psjsmodel.rate_variation:
                for codon_site_pair in self.data.space_codon_site_pair[n]:
                    ll, idx_pairs = self._sitewise_loglikelihood(n, codon_site_pair)
                    for i in range(len(ll)):
                        pair_to_lnL[idx_pairs[i]] = ll[i]
            else:
                ll, idx_pairs = self._sitewise_loglikelihood(n, codon_site_pair = None)
                for i in range(len(ll)):
                    pair_to_lnL[idx_pairs[i]] = ll[i]

        return pair_to_lnL

    def get_pairwise_loglikelihood_summary(self, summary_file):
        pair_to_lnL = self._sitewise_loglikelihood_for_all_n()
        with open(summary_file, 'w+') as f:
            f.write('#Site1\tSite2\tlnL\t\n')
            for pair in sorted(pair_to_lnL.keys()):
                f.write('\t'.join([str(pair[0]), str(pair[1]), str(pair_to_lnL[pair])]) + '\n')
            

    def loglikelihood_and_gradient_for_one_n(self, n):
        delta = 1e-8
        x = deepcopy(self.x)
        ll, edge_derivs = self._loglikelihood_for_one_n(n, edge_derivative = True)
            
        if self.root_by_dup:
            x_rate_derivs = edge_derivs
        else:
            two_rate_derivs = [[edge_derivs[i], self.tree.edge_list[i]] for i in range(len(edge_derivs)) if 'N0' in self.tree.edge_list[i]]
            assert(len(two_rate_derivs) == 2)  # should be exactly two branches
            merged_deriv = 0.0
            for i in range(2):
                if 'D1' in two_rate_derivs[i][1]:
                    merged_deriv += two_rate_derivs[i][0] * 10.0
                else:
                    merged_deriv += two_rate_derivs[i][0] * 10.0 / 9.0
                    
            x_rate_derivs = [merged_deriv] + \
                            [edge_derivs[i] for i in range(len(edge_derivs)) if not 'N0' in self.tree.edge_list[i]]

        m = len(self.x) - len(x_rate_derivs)

        # use finite differences to estimate derivatives with respect to these parameters
        other_derivs = []

        for i in range(m):
            if self.force == None or not i in self.force:
                x_plus_delta = np.array(self.x)
                x_plus_delta[i] += delta
                self.unpack_x(x_plus_delta)
                ll_delta, _ = self._loglikelihood_for_one_n(n, False)
                d_estimate = (ll_delta - ll) / delta
                other_derivs.append(d_estimate)
                # restore self.x
                self.unpack_x(x)
            else:
                other_derivs.append(0.0)
                
        other_derivs = np.array(other_derivs)
        f = -ll
        g = -np.concatenate((other_derivs, x_rate_derivs))

        return f, g
    
    def loglikelihood_and_gradient(self, display = False):
        sum_f = 0.0
        sum_g = 0.0
        inc = 0.05
        for it in range(len(self.data.space_list)):
            n = sorted(self.data.space_list)[it]
            if (it + 0.0) / len(self.data.space_list) > inc and display:
                print(str(floor(inc * 100)) + '% finished,  n =', n)
                inc += 0.2
            f, g = self.loglikelihood_and_gradient_for_one_n(n)
            sum_f += f
            sum_g += g

        # Now record lnL in self.ll
        self.ll = -f
        if display:
            print ('Sum log likelihood = ', sum_f)
            print ('All derivatives = ', sum_g)
            print ('Current exp x array = ', np.exp(self.x))
            print ('PM parameters = ' + ' '.join([i + ':' + str(self.psjsmodel.PMModel.parameters[i]) for i in self.psjsmodel.PMModel.parameter_list]))
            print ('Edge lengths = ', self.tree.edge_to_blen)
            print ()
            
        return sum_f, sum_g

    def objective_and_gradient(self,display, x):
        self.unpack_x(x)
        f, g = self.loglikelihood_and_gradient(display = display)
        self.auto_save += 1
        if self.auto_save == PSJSGeneconv.auto_save_step:
            self.save_x()
            self.auto_save = 0

        # Now save the lnL parameter and gradient values
        self.save_iteration(f, x, g)
        return f, g

    def save_iteration(self, f, x, df):
        with open(self.log_file, 'a') as g:
            g.write('\t'.join([str(item) for item in [f] + list(x) + list(df)]) + '\n')
            
    def objective_wo_gradient(self, display, x):
        self.unpack_x(x)

        sum_f = 0.0
        inc = 0.05
        for it in range(len(self.data.space_list)):
            n = sorted(self.data.space_list)[it]
            if (it + 0.0) / len(self.data.space_list) > inc and display:
                print(str(floor(inc * 100)) + '% finished,  n =', n)
                inc += 0.2
            f, g = self._loglikelihood_for_one_n(n, False)
            sum_f += f
            
        self.ll = sum_f
        if display:
            print ('log likelihood = ', sum_f)
            print ('Current x array = ', self.x)

        return -sum_f

    def objective_tract_p(self, display, log_tract_p):
        assert(self.psjsmodel.IGC_pm == 'One rate')
        # Only implemented for One rate model for now
        log_tau = self.psjsmodel.x_IGC[0] - self.psjsmodel.x_IGC[1]
        new_x_IGC = [log_tau + log_tract_p, log_tract_p]
        if self.root_by_dup:
            x_rate = self.x[-len(self.tree.edge_list):]
            x_js   = self.x[:(-len(self.tree.edge_list))]
        else:
            x_rate = self.x[-(len(self.tree.edge_list) - 1):]
            x_js   = self.x[:-(len(self.tree.edge_list) - 1)]

        x_js[-1] = new_x_IGC[-1]
        x_js[-2] = new_x_IGC[-2]
        return self.objective_wo_gradient(display, np.concatenate((x_js, x_rate)))

    def objective_2d_x_IGC(self, display, new_x_IGC):
        assert(self.psjsmodel.IGC_pm == 'One rate')
        # Only implemented for One rate model for now
        if self.root_by_dup:
            x_rate = self.x[-len(self.tree.edge_list):]
            x_js   = self.x[:(-len(self.tree.edge_list))]
        else:
            x_rate = self.x[-(len(self.tree.edge_list) - 1):]
            x_js   = self.x[:-(len(self.tree.edge_list) - 1)]

        x_js[-1] = new_x_IGC[-1]
        x_js[-2] = new_x_IGC[-2]
        return self.objective_wo_gradient(display, np.concatenate((x_js, x_rate)))

        
      

    def optimize_x_IGC(self, display = True, dimension = 1, method = 'L-BFGS-B'):
        if dimension == 1:
            f = partial(self.objective_tract_p, display)
            guess_x = self.psjsmodel.x_IGC[1]
            if method == 'L-BFGS-B':
                bnds = [(None, 0.0)]
            elif method == 'BasinHopping':
                bnds = [(-20.0, 0.0)]
            elif method == 'DifferentialEvolution':
                bnds = [(-20.0, 0.0)]
            else:
                sys.exit('Optimization method is not implemented!')
        if dimension == 2:
            f = partial(self.objective_2d_x_IGC, display)
            guess_x = self.psjsmodel.x_IGC
            if method == 'L-BFGS-B':
                bnds = [(None, None), (None, 0.0)]
            elif method == 'BasinHopping':
                bnds = [(-20.0, 20.0), (-20.0, 0.0)]
            elif method == 'DifferentialEvolution':
                bnds = [(-20.0, 20.0), (-20.0, 0.0)]
            else:
                sys.exit('Optimization method is not implemented!')
        if method == 'L-BFGS-B':
            result = scipy.optimize.minimize(f, guess_x, jac = False, method = method, bounds = bnds)
        elif method == 'BasinHopping':
            sys.exit('Not implemented yet!')
        else:
            sys.exit('Not implemented yet!')
            
        self.save_x()
        print (result)
        return result
            
    def get_mle(self, display = True, derivative = True):
        self.unpack_x(self.x)  # do one more update first
        if derivative:
            f = partial(self.objective_and_gradient, display)
        else:
            f = partial(self.objective_wo_gradient, display)

        guess_x = self.x
        bnds = [(None, -0.001)] * 3
        bnds.extend([(None, None)] * (len(self.x) - 3))
##        if self.root_by_dup:
##            bnds  = [(None, None)] * len(self.tree.edge_list)
##        else:
##            bnds  = [(None, None)] * (len(self.tree.edge_list) - 1)



        if derivative:
            result = scipy.optimize.minimize(f, guess_x, jac = True, method = 'L-BFGS-B', bounds = bnds)
        else:
            result = scipy.optimize.minimize(f, guess_x, jac = False, method = 'L-BFGS-B', bounds = bnds)

        self.save_x()
        print(result)
        return result

    def save_x(self):
        save = self.x
        np.savetxt(open(self.save_file, 'w+'), save.T)

    def initialize_by_save(self):
        self.x = np.loadtxt(open(self.save_file, 'r'))
        self.unpack_x(self.x)
        
    def get_summary(self):
        if self.ll == None:
            #self.ll = -1.0
            self.objective_wo_gradient(False, self.x)
            
        summary_mat = [self.ll, self.data.nsites]
        label = ['ll', 'length']
        
        for par in self.psjsmodel.PMModel.parameter_list:
            label.append(par)
            summary_mat.append(self.psjsmodel.PMModel.parameters[par])

        label.extend(['init_rate', 'tract_length'])
        summary_mat.extend([np.exp(self.psjsmodel.x_IGC[0]), 1.0 / np.exp(self.psjsmodel.x_IGC[1])])


        # TODO: Need to create a PSIGCModel class
##        for par in self.psjsmodel.IGCModel.parameter_list:
##            label.append(par)
##            summary_mat.append(self.psjsmodel.IGCModel.parameters[par])
        
        for edge in self.tree.edge_list:
            summary_mat.append(self.tree.edge_to_blen[edge])
            label.append('__'.join(edge))

        if not self.ExpectedGeneconv == None:
            for edge in self.tree.edge_list:
                label.append('__'.join(edge) + '__numIGC')
                summary_mat.append(self.ExpectedGeneconv[edge])

        return summary_mat, label

    def get_individual_summary(self, summary_file):
        summary, label = self.get_summary()
        summary = np.matrix(summary)
        footer = ' '.join(label)  # row labels

        np.savetxt(open(summary_file, 'w+'), summary.T, delimiter = ' ', footer = footer)
                    
        


if __name__ == '__main__':

#######################
    ###########Yeast
#######################
    gene_to_orlg_file = '../test/YDR418W_YEL054C_GeneToOrlg.txt'
    alignment_file = '../test/YDR418W_YEL054C_test.fasta'

    tree_newick = '../test/YeastTree.newick'
    DupLosList = '../test/YeastTestDupLost.txt'
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    seq_index_file = '../test/YDR418W_YEL054C_seq_index.txt'
    seq_index_file = None

    pm_model = 'HKY'
    
    IGC_pm = 'One rate'
    space_list = range(1, 330, 20)
    space_list = None

    
    cdna = False
    allow_same_codon = False
    rate_variation = False
    alignment_file = '../test/YDR418W_YEL054C_input.fasta'
    save_file = '../test/save/PSJS_HKY_YDR418W_YEL054C_nonclock_save.txt'
    log_file  = '../test/log/PSJS_HKY_YDR418W_YEL054C_nonclock_log.txt'
    summary_file = '../test/Summary/PSJS_HKY_YDR418W_YEL054C_nonclock_summary.txt'
    x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3, 1.0 / 30.0 ])
    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, cdna, allow_same_codon, tree_newick, DupLosList,x_js, pm_model, IGC_pm, rate_variation,
                      node_to_pos, terminal_node_list, save_file, log_file, space_list = space_list)
##    scene = test.get_scene(469, None)

##    alignment_file = '../test/YDR418W_YEL054C_input.fasta'
##    cdna = True
##    allow_same_codon = True
##    rate_variation = True
##    save_file = '../test/save/PSJS_HKY_YDR418W_YEL054C_rv_nonclock_save.txt'
##    log_file  = '../test/log/PSJS_HKY_YDR418W_YEL054C_rv_nonclock_log.txt'
##    summary_file = '../test/Summary/PSJS_HKY_YDR418W_YEL054C_rv_nonclock_summary.txt'
##    x_js = np.log([ 0.4, 0.6, 0.7,  4.35588244, 0.8, 1.8,  0.3, 1.0 / 30.0 ])
##    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, cdna, allow_same_codon, tree_newick, DupLosList,x_js, pm_model, IGC_pm, rate_variation,
##                      node_to_pos, terminal_node_list, save_file, log_file, space_list = space_list)
##
##    x = np.array([ -0.92969299,  -0.60831311,  -1.00401104,  -0.8891694 ,
##        -0.89092996,  -0.20694348, -65.92436066, 0.0, -2.73221044,
##       -26.7831379 ,  -3.44536047,  -2.72202154, -21.09784191,
##        -2.70441861, -22.63017035,  -2.72779708, -20.53753293,
##        -3.44335198,  -2.73358356, -20.14414384])
##    test.unpack_x(x)
    
##    self = test
##    test.cal_iid_observations()
##    print (test.iid_observations.keys() )
##    scene1 = test.get_scene(13, (1, 2))
##    scene2 = test.get_scene(13, (2, 3))
##    scene3 = test.get_scene(13, (3, 1))
##    for key in scene1:
##        try:
##            print (key, scene1[key] == scene2[key], scene1[key] == scene3[key])
##        except:
##            print (key)
##    key = 'process_definitions'
##    print (key, scene1[key][0] == scene2[key][0], scene1[key][0] == scene3[key][0])
##    print (key, scene1[key][1] == scene2[key][1], scene1[key][1] == scene3[key][1])
##    print(test._loglikelihood(13,  codon_site_pair = (2, 3), edge_derivative = False), test._loglikelihood(13,  codon_site_pair = (1, 2), edge_derivative = False),
##          test._loglikelihood(13,  codon_site_pair = (3, 1), edge_derivative = False))
##    print(sum([test._loglikelihood(13,  codon_site_pair = (2, 3), edge_derivative = False)[0], test._loglikelihood(13,  codon_site_pair = (1, 2), edge_derivative = False)[0],
##          test._loglikelihood(13,  codon_site_pair = (3, 1), edge_derivative = False)[0]]))

##    print(test._loglikelihood_for_one_n(13, False))
    #print(test.objective_wo_gradient(True, test.x))
    #test.optimize_x_IGC(True)
    #print(test.loglikelihood_and_gradient_for_one_n(13))
    pairwise_lnL_summary_file = summary_file.replace('_summary.txt', '_lnL_summary.txt')
    test.get_pairwise_loglikelihood_summary(pairwise_lnL_summary_file)
##    ll, g = test.loglikelihood_and_gradient(True)
##    print( ll/(test.data.nsites - 1.0), g/(test.data.nsites - 1.0) )
##    for edge in test.tree.edge_list:
##        print (edge, test.tree.edge_to_blen[edge])
    #print(test.loglikelihood_and_gradient_for_one_n(n))
    #test.get_mle()
    #test.get_individual_summary(summary_file)


#########################
##    ###########EDN ECP
#########################
##    gene_to_orlg_file = '../test/EDN_ECP_GeneToOrlg.txt'
##    alignment_file = '../test/EDN_ECP_Cleaned.fasta'
##    tree_newick = '../test/input_tree.newick'
##    DupLosList = '../test/EDN_ECPDupLost.txt'
##    terminal_node_list = ['Chimpanzee', 'Gorilla', 'Orangutan', 'Macaque', 'Tamarin']
##    node_to_pos = {'D1':0}
##    seq_index_file = None
##    pm_model = 'HKY'
##    IGC_pm = 'One rate'
##    space_list = None
##
##    cdna = True
##    allow_same_codon = True
##    rate_variation = False
##    save_file = '../test/save/PSJS_HKY_EDN_ECP_nonclock_save.txt'
##    summary_file = '../test/Summary/PSJS_HKY_EDN_ECP_nonclock_summary.txt'
##    x_js = np.log([ 0.4, 0.6, 0.7,  4.35588244,  0.3, 1.0 / 30.0 ])
##    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, cdna, allow_same_codon, tree_newick, DupLosList,x_js, pm_model, IGC_pm, rate_variation,
##                      node_to_pos, terminal_node_list, save_file, space_list = space_list)
##    #test.get_mle()
##    #test.get_individual_summary(summary_file)
##    pairwise_lnL_summary_file = '../test/Summary/PSJS_HKY_EDN_ECP_nonclock_lnL_summary.txt'
##    test.get_pairwise_loglikelihood_summary(pairwise_lnL_summary_file)
##    #ll, idx_pairs = test._sitewise_loglikelihood(467, (1, 3))
    
    






