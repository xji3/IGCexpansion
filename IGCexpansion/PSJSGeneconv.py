# A separate file for JSLikelihood
# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Xiang Ji
# xji3@ncsu.edu
from __future__ import print_function
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
                 
                 tree_newick, DupLosList,                 # Tree input
                 x_js, pm_model, IGC_pm,                  # JSModel input
                 node_to_pos, terminal_node_list,         # Configuration input
                 save_file,                               # Auto save file
#                root_by_dup = False,                     # JSModel input
                 force = None,                            # Parameter value constraint
                 nsites = None, space_list = None
                 ):  

        
        self.tree = Tree(tree_newick, DupLosList, terminal_node_list, node_to_pos)
        self.data = Data(alignment_file, gene_to_orlg_file, two_sites = True, space_list = space_list)
        
        self.psjsmodel = PSJSModel(x_js, pm_model, self.tree.n_orlg, IGC_pm, force)
        self.node_to_pos = node_to_pos
        self.terminal_node_list = terminal_node_list
        self.root_by_dup = self.tree.is_duplication_node(self.tree.phylo_tree.root.name)
        self.save_file = save_file
        self.force = force
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
    
    def get_scene(self, n):
        state_space_shape = self.psjsmodel.state_space_shape
        conf_list = count_process(self.tree.node_to_conf)
        process_definitions = [self.psjsmodel.get_PM_process_definition(), self.psjsmodel.get_IGC_process_definition(n)]
        self.tree.get_tree_process(conf_list)

        prior_feasible_states, distn = self.get_prior()
        
        if self.nsites is None:
            observable_nodes, observable_axes, iid_observations = get_PS_iid_observations(self.data, self.tree, self.data.nsites - n, n, self.psjsmodel.PMModel.data_type)
        else:
            observable_nodes, observable_axes, iid_observations = get_PS_iid_observations(self.data, self.tree, self.nsites, n, self.psjsmodel.PMModel.data_type)
        scene = dict(
            node_count = len(self.tree.node_to_num),
            process_count = len(process_definitions),
            state_space_shape = state_space_shape,
            tree = self.tree.tree_json,
            root_prior = {'states':prior_feasible_states,
                          'probabilities':distn},
            process_definitions = process_definitions,
            observed_data = {
                'nodes':observable_nodes,
                'variables':observable_axes,
                'iid_observations':iid_observations
                }
            )
        return scene


    def _loglikelihood(self, n, edge_derivative = False):
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
        self.ll = ll
        if edge_derivative:
            edge_derivs = j_out['responses'][1]
        else:
            edge_derivs = []

        return ll, edge_derivs

    def loglikelihood_and_gradient_for_one_n(self, n):
        delta = 1e-8
        x = deepcopy(self.x)
        ll, edge_derivs = self._loglikelihood(n, edge_derivative = True)
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
                ll_delta, _ = self._loglikelihood(n, False)
                d_estimate = (ll_delta - ll) / delta
                other_derivs.append(d_estimate)
                # restore self.x
                self.unpack_x(x)
            else:
                other_derivs.append(0.0)
                
        other_derivs = np.array(other_derivs)
        self.ll = ll
        f = -ll
        g = -np.concatenate((other_derivs, x_rate_derivs))
        if False:
            print ('log likelihood = ', ll)
            print ('Edge derivatives = ', x_rate_derivs)
            print ('other derivatives:', other_derivs)
            print ('Current exp x array = ', np.exp(self.x))
            print ('PM parameters = ' + ' '.join([i + ':' + str(self.psjsmodel.PMModel.parameters[i]) for i in self.psjsmodel.PMModel.parameter_list]))
            print ('Edge lengths = ', self.tree.edge_to_blen)
            print ()
            
        return f, g
    
    def loglikelihood_and_gradient(self, display = False):
        sum_f = 0.0
        sum_g = 0.0
        inc = 0.05
        for it in range(len(self.data.two_sites_name_to_seq.keys())):
            n = sorted(self.data.two_sites_name_to_seq.keys())[it]
            if (it + 0.0) / len(self.data.two_sites_name_to_seq) > inc and display:
                print(str(floor(inc * 100)) + '% finished,  n =', n)
                inc += 0.2
            f, g = self.loglikelihood_and_gradient_for_one_n(n)
            sum_f += f
            sum_g += g
            
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
        return f, g

    def objective_wo_gradient(self, display, x):
        self.unpack_x(x)
        sum_f = 0.0
        for n in self.data.two_sites_name_to_seq.keys():
            f, g = self._loglikelihood(n, False)
            sum_f += f

        if display:
            print ('log likelihood = ', sum_f)
            print ('Current x array = ', self.x)

        return -sum_f
            
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
        summary_mat = [self.ll]
        label = ['ll']
        
        for par in self.jsmodel.PMModel.parameter_list:
            label.append(par)
            summary_mat.append(self.jsmodel.PMModel.parameters[par])

        for par in self.jsmodel.IGCModel.parameter_list:
            label.append(par)
            summary_mat.append(self.jsmodel.IGCModel.parameters[par])
        
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
    alignment_file = '../test/YDR418W_YEL054C_MG94_geo_10.0_Sim_8.fasta'

    tree_newick = '../test/YeastTree.newick'
    DupLosList = '../test/YeastTestDupLost.txt'
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    save_file = '../test/save/PSJS_HKY_YDR418W_YEL054C_nonclock_save.txt'
    summary_file = '../test/Summary/PSJS_HKY_YDR418W_YEL054C_nonclock_summary.txt'

    pm_model = 'HKY'
    x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3, 1.0 / 30.0 ])
    IGC_pm = 'One rate'
    space_list = range(1, 330, 20)
    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList,x_js, pm_model, IGC_pm,
                      node_to_pos, terminal_node_list, save_file, space_list = space_list)
    self = test
    #print(test.loglikelihood_and_gradient_for_one_n(n))
    #test.get_mle()







