# A separate file for JSLikelihood
# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Xiang Ji
# xji3@ncsu.edu
from __future__ import print_function
import jsonctmctree.ll, jsonctmctree.interface
from Data import Data
from Tree import Tree
from JSModel import JSModel
from Func import *
from copy import deepcopy
from functools import partial
import scipy
import scipy.optimize
import os
from Common import *
from math import floor

from TriGeneconv import *

class JSGeneconv:
    auto_save_step = 2
    def __init__(self, alignment_file, gene_to_orlg_file, cdna, # Data input
                 tree_newick, DupLosList,                       # Tree input
                 x_js, pm_model, IGC_pm, rate_variation,        # JSModel input
                 node_to_pos, terminal_node_list,               # Configuration input
                 save_file,                                     # Auto save file
                 force = None,                                  # Parameter value constraint
                 seq_index_file = None,                         # Seq index file
                 nsites = None):  

        
        self.tree = Tree(tree_newick, DupLosList, terminal_node_list, node_to_pos)
        self.data = Data(alignment_file, gene_to_orlg_file, seq_index_file = seq_index_file, cdna = cdna)
        # added new variable for jsmodel
        conf_list = count_process(self.tree.node_to_conf)
        accessible_orlg_pair = get_accessible_orlg_pair(conf_list)
        
        self.jsmodel = JSModel(self.tree.n_js, x_js, pm_model, self.tree.n_orlg, IGC_pm, rate_variation, accessible_orlg_pair, force)
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
            self.unpack_x(self.x)

        assert(self.jsmodel.n_orlg == self.tree.n_orlg)  # make sure n_orlg got updated
        self.ll   = None              # used to store log likelihood
        self.ExpectedGeneconv = None  # Used to store expectedNumGeneconv
        self.auto_save = 0            # used to control auto save frequency

        self.nsites = nsites
        self.iid_observations = None  # Store iid_observations which only needs one time calculation
        self.observable_nodes = None
        self.observable_axes  = None

        assert(self.self_check())
        

            
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
        self.jsmodel.update_by_x_js(x_js)

    def self_check(self):
        check_status = True
        if self.jsmodel.rate_variation:
            check_status = check_status and self.data.cdna

        return check_status

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
        configuration = deepcopy(self.tree.node_to_conf[self.tree.phylo_tree.root.name])
        # assign prior feasible states and its distribution with given configuration
        if self.tree.phylo_tree.root.name in self.tree.node_to_dup:
            assert(self.root_by_dup)
            new_orlgs = self.tree.node_to_dup[self.tree.phylo_tree.root.name]
            for i in range(len(configuration)):
                if configuration[i][0] in new_orlgs:
                    configuration[i][0] = new_orlgs[0]
            
        ortho_group_to_pos = divide_configuration(configuration)
        assert( not ortho_group_to_pos['distinct']) # don't allow distinct lineages at root for now
        assert( len(ortho_group_to_pos['extent']) < 3) # root is on either outgroup lineage or is the first duplication node
        
        prior_feasible_states = []
        distn = []
        for nt in range(self.jsmodel.state_space_shape[0]):
            js_state = [nt] * self.jsmodel.n_js
            if self.jsmodel.is_state_compatible(js_state, configuration):
                prior_feasible_states.append(js_state)
                distn.append(self.jsmodel.PMModel.get_stationary_distn(nt))

        distn = np.array(distn) / sum(distn)
        return prior_feasible_states, distn
    
    def get_scene(self):
        state_space_shape = self.jsmodel.state_space_shape
        conf_list = count_process(self.tree.node_to_conf)
        
        self.tree.get_tree_process(conf_list)

        prior_feasible_states, distn = self.get_prior()

        if self.iid_observations == None:
            self.cal_iid_observations()
            
        if self.jsmodel.rate_variation:
            scene_list = []
            for codon_site in range(1, 4):
                process_definitions, conf_list = get_process_definitions(self.tree, self.jsmodel, codon_site = codon_site)
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
                        'iid_observations':self.iid_observations[codon_site]
                        } 
                    )
                scene_list.append(scene)
            return scene_list
        else:
            process_definitions, conf_list = get_process_definitions(self.tree, self.jsmodel)
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
                    'iid_observations':self.iid_observations
                    }
                )
            return scene

    def cal_iid_observations(self):
        if self.iid_observations == None:
            if self.data.cdna:
                self.observable_nodes, self.observable_axes, self.iid_observations = get_iid_observations(self.data, self.tree, self.nsites, self.jsmodel.PMModel.data_type)
                if not self.jsmodel.rate_variation:
                    self.iid_observations = [i for j in range(1,4) for i in self.iid_observations[j]]
            else:
                if self.nsites is None:
                    self.observable_nodes, self.observable_axes, self.iid_observations = get_iid_observations(self.data, self.tree, self.data.nsites, self.jsmodel.PMModel.data_type)
                else:
                    self.observable_nodes, self.observable_axes, self.iid_observations = get_iid_observations(self.data, self.tree, self.nsites, self.jsmodel.PMModel.data_type)


        
    def get_expectedNumGeneconv(self, display = False):
        process_definitions, conf_list = get_process_definitions(self.tree, self.jsmodel, proportions = True)
        requests = [{'property' : 'SDNTRAN', 'transition_reduction' : process_definitions[i]} for i in range(len(process_definitions))]
        scene = self.get_scene()
        j_in = {'scene' : scene,
                'requests':requests}
        j_out = jsonctmctree.interface.process_json_in(j_in)

        status = j_out['status']
        self.ExpectedGeneconv = {self.tree.edge_list[i] : j_out['responses'][scene['tree']['edge_processes'][i]][i] for i in range(len(self.tree.edge_list))}


        
    def get_pairDirectionalExpectedNumGeneconv(self, orlg_pair, display = False):
        process_definitions, conf_list = get_directional_process_definitions(self.tree, self.jsmodel, orlg_pair)
        requests = [{'property' : 'SDNTRAN', 'transition_reduction' : process_definitions[i]} for i in range(len(process_definitions))]
        scene = self.get_scene()
        j_in = {'scene' : scene,
                'requests':requests}
        j_out = jsonctmctree.interface.process_json_in(j_in)

        status = j_out['status']
        pairDirectionalExpectedGeneconv = {self.tree.edge_list[i] : j_out['responses'][scene['tree']['edge_processes'][i]][i] for i in range(len(self.tree.edge_list))}
        return pairDirectionalExpectedGeneconv 



    def _loglikelihood(self, edge_derivative = False):
        scene = self.get_scene()
        log_likelihood_request = {'property':'snnlogl'}
        derivatives_request = {'property':'sdnderi'}
        if edge_derivative:
            requests = [log_likelihood_request, derivatives_request]
        else:
            requests = [log_likelihood_request]

        if self.jsmodel.rate_variation:
            edge_derivs = 0.0
            ll = 0.0
            for codon_site_scene in scene:
                j_in = {
                    'scene' : codon_site_scene,
                    'requests' : requests
                    }
                j_out = jsonctmctree.interface.process_json_in(j_in)
                status = j_out['status']
                ll += j_out['responses'][0]

                if edge_derivative:
                    #print('In _loglikelihood() edge_derivs: ', j_out['responses'][1])
                    edge_derivs += np.array(j_out['responses'][1])
                else:
                    edge_derivs = []
        else:
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

    def _sitewise_loglikelihood(self):
        scene = self.get_scene()
        log_likelihood_request = {'property':'dnnlogl'}
        requests = [log_likelihood_request]


        if self.jsmodel.rate_variation:
            edge_derivs = 0.0
            ll_list = []
            for codon_site_scene in scene:
                j_in = {
                    'scene' : codon_site_scene,
                    'requests' : requests
                    }
                j_out = jsonctmctree.interface.process_json_in(j_in)
                status = j_out['status']
                ll_list.append(j_out['responses'][0])
            ll = []
            for i in range(self.data.nsites):
                codon_site = i % 3
                pos = int(floor((i + 0.5) / 3))
                ll.append(ll_list[codon_site][pos])
        
        else:
            j_in = {
                'scene' : scene,
                'requests' : requests
                }
            j_out = jsonctmctree.interface.process_json_in(j_in)

            status = j_out['status']
        
            ll = j_out['responses'][0]

        return ll

    def get_sitewise_loglikelihood_summary(self, summary_file):
        ll = self._sitewise_loglikelihood()
        with open(summary_file, 'w+') as f:
            f.write('#Site\tlnL\t\n')
            for i in range(self.data.nsites):
                f.write('\t'.join([str(i), str(ll[i])]) + '\n')
       

    
    def loglikelihood_and_gradient(self, display = False):
        delta = 1e-8
        x = deepcopy(self.x)
        ll, edge_derivs = self._loglikelihood(edge_derivative = True)
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
                x_plus_delta = deepcopy(np.array(self.x))
                x_plus_delta[i] += delta
                self.unpack_x(x_plus_delta)
                ll_plus, _ = self._loglikelihood(False)
                x_plus_delta[i] -= 2*delta
                self.unpack_x(x_plus_delta)
                ll_minus, _ = self._loglikelihood(False)
                d_estimate = (ll_plus - ll_minus) / (2 * delta)
                other_derivs.append(d_estimate)
                # restore x_plut_delta
                x_plus_delta[i] += delta
                # restore self.x
                self.unpack_x(x)
            else:
                other_derivs.append(0.0)
                
        other_derivs = np.array(other_derivs)
        self.ll = ll
        f = -ll
        g = -np.concatenate((other_derivs, x_rate_derivs))
        if display:
            print ('log likelihood = ', ll)
            print ('Edge derivatives = ', x_rate_derivs)
            print ('other derivatives:', other_derivs)
            print ('Current exp x array = ', np.exp(self.x))
            print ('PM parameters = ' + ' '.join([i + ':' + str(self.jsmodel.PMModel.parameters[i]) for i in self.jsmodel.PMModel.parameter_list]))
            print ('IGC parameters = ' + ' '.join([i + ':' + str(self.jsmodel.IGCModel.parameters[i]) for i in self.jsmodel.IGCModel.parameter_list]))
            print ('Edge lengths = ', self.tree.edge_to_blen)
            print ()
            
        return f, g

    def objective_and_gradient(self,display, x):
        self.unpack_x(x)
        f, g = self.loglikelihood_and_gradient(display = display)
        self.auto_save += 1
        if self.auto_save == JSGeneconv.auto_save_step:
            self.save_x()
            self.auto_save = 0
        return f, g

    def objective_wo_gradient(self, display, x):
        self.unpack_x(x)
        ll = self._loglikelihood()[0]
        self.ll = ll
        if display:
            print ('log likelihood = ', ll)
            print ('Current x array = ', self.x)

        return -ll
            
    def get_mle(self, display = True, derivative = True, method = 'L-BFGS-B', niter = 2000):
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


        if method == 'L-BFGS-B':
            if derivative:
                result = scipy.optimize.minimize(f, guess_x, jac = True, method = 'L-BFGS-B', bounds = bnds)
            else:
                result = scipy.optimize.minimize(f, guess_x, jac = False, method = 'L-BFGS-B', bounds = bnds)
        elif method == 'BasinHopping':
            bnds = bnds = [(-20.0, -0.001)] * 3
            bnds.extend([(-20.0, 20.0)] * 2 + [(-20.0, 2.0)] * (len(self.x) - 5))
            if derivative:
                result = scipy.optimize.basinhopping(f, guess_x, minimizer_kwargs = {'method':'L-BFGS-B', 'jac':True, 'bounds':bnds}, niter = niter)#, callback = self.check_boundary)
            else:
                result = scipy.optimize.basinhopping(f, guess_x, minimizer_kwargs = {'method':'L-BFGS-B', 'jac':False, 'bounds':bnds}, niter = niter)#, callback = self.check_boundary)
        elif method == 'DifferentialEvolution':
            f = partial(self.objective_wo_gradient, display)
            bnds = bnds = [(-20.0, -0.001)] * 3
            bnds.extend([(-20.0, 20.0)] * 2 + [(-20.0, 2.0)] * (len(self.x) - 5))
            result = scipy.optimize.differential_evolution(f, bnds)
            

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
        summary_mat = [self.ll, self.data.nsites]
        label = ['ll', 'length']
        
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
    ###########EDN ECP
#######################
    gene_to_orlg_file = '../test/EDN_ECP_GeneToOrlg.txt'
    alignment_file = '../test/EDN_ECP_Cleaned.fasta'
    tree_newick = '../test/input_tree.newick'
    DupLosList = '../test/EDN_ECPDupLost.txt'
    terminal_node_list = ['Chimpanzee', 'Gorilla', 'Orangutan', 'Macaque', 'Tamarin']
    node_to_pos = {'D1':0}
    seq_index_file = None
    pm_model = 'HKY'
    IGC_pm = 'One rate'
    space_list = None

    cdna = True
    allow_same_codon = True
    rate_variation = False
    save_file = '../test/save/JS_HKY_EDN_ECP_nonclock_save.txt'
    summary_file = '../test/Summary/JS_HKY_EDN_ECP_nonclock_summary.txt'
    x_js = np.log([ 0.4, 0.6, 0.7,  4.35588244,  0.3])
    force = None
    test = JSGeneconv(alignment_file, gene_to_orlg_file, cdna, tree_newick, DupLosList,x_js, pm_model, IGC_pm, rate_variation,
                      node_to_pos, terminal_node_list, save_file, force)
    test.get_mle()
    test.get_individual_summary(summary_file)


    cdna = True
    allow_same_codon = True
    rate_variation = False
    save_file = '../test/save/Force_JS_HKY_EDN_ECP_nonclock_save.txt'
    summary_file = '../test/Summary/Force_JS_HKY_EDN_ECP_nonclock_summary.txt'
    x_js = np.log([ 0.4, 0.6, 0.7,  4.35588244,  0.3])
    force = {4:0.0}
    test = JSGeneconv(alignment_file, gene_to_orlg_file, cdna, tree_newick, DupLosList,x_js, pm_model, IGC_pm, rate_variation,
                      node_to_pos, terminal_node_list, save_file, force)
    test.get_mle()
    test.get_individual_summary(summary_file)

    cdna = True
    allow_same_codon = True
    rate_variation = True
    save_file = '../test/save/JS_HKY_rv_EDN_ECP_nonclock_save.txt'
    summary_file = '../test/Summary/JS_HKY_rv_EDN_ECP_nonclock_summary.txt'
    x_js = np.log([ 0.4, 0.6, 0.7,  4.35588244, 0.8, 9.0,  0.3])
    force = None
    test = JSGeneconv(alignment_file, gene_to_orlg_file, cdna, tree_newick, DupLosList,x_js, pm_model, IGC_pm, rate_variation,
                      node_to_pos, terminal_node_list, save_file, force)
    test.get_mle()
    test.get_individual_summary(summary_file)


    cdna = True
    allow_same_codon = True
    rate_variation = True
    save_file = '../test/save/Force_JS_HKY_rv_EDN_ECP_nonclock_save.txt'
    summary_file = '../test/Summary/Force_JS_HKY_rv_EDN_ECP_nonclock_summary.txt'
    x_js = np.log([ 0.4, 0.6, 0.7,  4.35588244, 0.8, 9.0,  0.3])
    force = {6:0.0}
    test = JSGeneconv(alignment_file, gene_to_orlg_file, cdna, tree_newick, DupLosList,x_js, pm_model, IGC_pm, rate_variation,
                      node_to_pos, terminal_node_list, save_file, force)
    test.get_mle()
    test.get_individual_summary(summary_file)
#########################
##    ###########Yeast
#########################
##    gene_to_orlg_file = '../test/YDR418W_YEL054C_GeneToOrlg.txt'
##    alignment_file = '../test/YDR418W_YEL054C_input.fasta'
##    cdna = True
##
##    tree_newick = '../test/YeastTree.newick'
##    DupLosList = '../test/YeastTestDupLost.txt'
##    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
##    node_to_pos = {'D1':0}
##    pm_model = 'HKY'
##    IGC_pm = 'One rate'
##
##
##
####    save_file = '../test/save/HKY_YDR418W_YEL054C_nonclock_One_rate_save.txt'
####    summary_file = '../test/Summary/HKY_YDR418W_YEL054C_nonclock_One_rate_summary.txt'    
####    rate_variation = False
####    x_js = np.log([ 0.1,   0.7,   0.1,  4.35588244,   1.01054376])    
####    test_old = JSGeneconv(alignment_file, gene_to_orlg_file, cdna, tree_newick, DupLosList,x_js, pm_model, IGC_pm, rate_variation,
####                      node_to_pos, terminal_node_list, save_file)
####    test_old.get_mle()
####    print (test_old._loglikelihood(True))
####
####    save_file = '../test/save/Force_HKY_YDR418W_YEL054C_nonclock_One_rate_save.txt'
####    summary_file = '../test/Summary/Force_HKY_YDR418W_YEL054C_nonclock_One_rate_summary.txt'    
####    rate_variation = False
####    force = {4:0.0}
####    x_js = np.log([ 0.1,   0.7,   0.1,  4.35588244,   1.01054376])    
####    test_old_force = JSGeneconv(alignment_file, gene_to_orlg_file, cdna, tree_newick, DupLosList,x_js, pm_model, IGC_pm, rate_variation,
####                      node_to_pos, terminal_node_list, save_file, force)
####    test_old_force.get_mle()
####    print (test_old_force._loglikelihood(True))
##
##
##    alignment_file = '../test/YDR418W_YEL054C_test.fasta'
##    alignment_file = '../test/YDR418W_YEL054C_input.fasta'
##    save_file = '../test/save/HKY_YDR418W_YEL054C_nonclock_rv_One_rate__save.txt'
##    summary_file = '../test/Summary/HKY_YDR418W_YEL054C_nonclock_rv_One_rate_summary.txt'
##    rate_variation = True
##    x_js = np.log([0.3, 0.5, 0.2, 9.5, 1.0, 2.6, 5.9])
##    #x_js = np.concatenate((test_old.jsmodel.x_js[:-1], [0.0, 0.0, test_old.jsmodel.x_js[-1]]))
##    #x = np.concatenate((x_js, test_old.x[len(test_old.jsmodel.x_js):]))
##    test = JSGeneconv(alignment_file, gene_to_orlg_file, cdna, tree_newick, DupLosList,x_js, pm_model, IGC_pm, rate_variation,
##                      node_to_pos, terminal_node_list, save_file)
##    #test.get_mle()
##    x = np.array([ -0.92969299,  -0.60831311,  -1.00401104,  -0.8891694 ,
##        -0.89092996,  -0.20694348, -65.92436066,  -2.73221044,
##       -26.7831379 ,  -3.44536047,  -2.72202154, -21.09784191,
##        -2.70441861, -22.63017035,  -2.72779708, -20.53753293,
##        -3.44335198,  -2.73358356, -20.14414384])
##    test.unpack_x(x)
##    print(test.x)
##    print(test.loglikelihood_and_gradient())
##    for edge in test.tree.edge_list:
##        print (edge, test.tree.edge_to_blen[edge])    
##
####    save_file = '../test/save/Force_HKY_YDR418W_YEL054C_nonclock_rv_One_rate__save.txt'
####    summary_file = '../test/Summary/Force_HKY_YDR418W_YEL054C_nonclock_rv_One_rate_summary.txt'
####
####    rate_variation = True
####    force = {6:0.0}
####    #x_js = np.concatenate((test_old.jsmodel.x_js[:-1], [0.0, 0.0, test_old.jsmodel.x_js[-1]]))
####    x_js = np.log([0.3, 0.5, 0.2, 9.5, 1.0, 2.6, 5.9])
####    test_force = JSGeneconv(alignment_file, gene_to_orlg_file, cdna, tree_newick, DupLosList,x_js, pm_model, IGC_pm, rate_variation,
####                      node_to_pos, terminal_node_list, save_file, force)
####    x = np.concatenate((x_js, test_old.x[len(test_old.jsmodel.x_js):]))
####    test_force.get_mle()
####
####    
####    self = test
####    print (test._loglikelihood(True))
####    test.cal_iid_observations()
####    scene = test.get_scene()
####    print(test.loglikelihood_and_gradient(True))
##
##    #test.get_individual_summary(summary_file)
##
##
##    #test.get_mle(method = 'BasinHopping')
####    test.get_mle(method = 'DifferentialEvolution')
####    self = test
####    print(test.tree.n_js, test.tree.n_orlg)
####
####    pm_model = 'HKY'
####    x_js = np.log([ 0.1,   0.7,   0.1,  4.35588244,   1.01054376, 1.5456])
####    IGC_pm = 'Most general'
####    save_file = '../test/save/HKY_YDR418W_YEL054C_nonclock_Most_general_save.txt'
####    summary_file = '../test/HKY_YLR406C_YDL075W_nonclock_dir_summary.txt'
####
####    test_2 = JSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList,x_js, pm_model, IGC_pm,
####                      node_to_pos, terminal_node_list, save_file)
####
####    self = test_2
###    print test.get_scene(100)
##    #test_2.get_mle(method = 'BasinHopping')
##    #test.get_expectedNumGeneconv()
##    #print(test.get_pairDirectionalExpectedNumGeneconv([1, 2]))
##    #test.get_individual_summary(summary_file)
##


###########################
####    ###########Class 1 ADH Genes
###########################
##    gene_to_orlg_file = '../test/Trigeneconv_GeneToOrlg.txt'
##    alignment_file = '../test/Trigeneconv_ADH_intron_input.fasta'
##    
##    data = Data(alignment_file, gene_to_orlg_file)
##    
##    tree_newick = '../test/Trigeneconv_ADH1Class_tree.newick'
##    DupLosList = '../test/Trigeneconv_ADH_DupLost.txt'
##    save_file = '../test/save/Trigeneconv_ADH_HKY_save.txt'
##
##    node_to_pos = {'D1':0, 'D2':0}
##    terminal_node_list = ['Baboon', 'Orangutan', 'Gorilla', 'Bonobo', 'Chimpanzee', 'Human']
##    tree = Tree(tree_newick, DupLosList, terminal_node_list, node_to_pos)
##    conf_list = count_process(tree.node_to_conf)
##    accessible_orlg_pair = get_accessible_orlg_pair(conf_list)
##    
##    pm_model = 'HKY'
##    x_pm = np.log([ 0.49978115,   0.60208772,   0.41240341,  5.35588244])
##    x_js = np.concatenate((x_pm, np.log(range(2, 2 + len(accessible_orlg_pair) * 2))))
##    IGC_pm = 'Most general'
##    nsites = None
##
##    save_file = '../test/save/Trigeneconv_ADH_HKY_dir_save.txt'
##
##    x_js = np.concatenate((x_pm, np.log(range(2, 2 + len(accessible_orlg_pair)))))
##    save_file = '../test/save/Trigeneconv_ADH_HKY_sym_save.txt'
##    IGC_pm = 'Symmetric general'
##    #force = {4:0.0}
##    force = None
##
##    test = JSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
##                      node_to_pos, terminal_node_list, save_file, nsites = nsites, force = force)
##    self = test
##    test.get_mle()
##    x_process = [-0.67216797,  -0.40212794,  -1.04889601,   1.20207224, -0.60249277]
##    x_rates = np.log([0.00030354302771, # D1_D2, N0_N1
##                      0.0613673848999,  # D2_N0, N1_N2
##                      0.010962330683,   # N0_N1, N2_N3
##                      0.0348850389125,  # N0_Ba, N2_Ba
##                      0.0131843207937,  # N1_Or, N3_Or
##                      0.00784978378118, # N1_N2, N3_N4
##                      0.00297723023713, # N2_N3, N4_N5
##                      0.00496736099364, # N2_Go, N4_Go
##                      0.000933915293812,# N3_Bo, N5_Bo
##                      8.05308221277e-07,# N3_N4, N5_N6
##                      0.00656242294756, # N4_Hu, N6_Hu
##                      0.00105781666618  # N4_Ch, N6_Ch
##                      ])
##    x_tri = np.concatenate((x_process, x_rates))
##    test.unpack_x(x_tri)
##    x = np.concatenate((x_js, np.log([0.2] * (len(test.tree.edge_list) ))))
##    test.unpack_x(x)
    #process_definitions, conf_list = get_process_definitions(self.tree, self.jsmodel)
    
##    conf_list = count_process(test.tree.node_to_conf)
##    configuration = conf_list[0]
###    process = jsmodel.get_process_definition(configuration)
##
##    edge_derivative = False
##    data = test.data
##    tree = test.tree
##    data_type = test.jsmodel.PMModel.data_type
##    
##    observable_nodes, observable_axes, iid_observations = get_iid_observations(test.data, test.tree, test.data.nsites, test.jsmodel.PMModel.data_type)
##    print (test._loglikelihood())
##    ExpectedGeneconv = test.get_expectedNumGeneconv()
##    pairExGeneconv = test.get_pairDirectionalExpectedNumGeneconv([3, 2])
##    reverse_pairExGeneconv = test.get_pairDirectionalExpectedNumGeneconv([2, 3])
##    print (pairExGeneconv)
##    print (reverse_pairExGeneconv)
##    #test.get_mle()
##
##
########from TriGeneconv.py
##    alignment_file = '../../TriGeneconv/data/ADH_intron_input.fasta'
##    newicktree = '../../TriGeneconv/data/ADH1Class_tree.newick'
##    save_path = '../../TriGeneconv/save/'
##    paralog = ['ADH1A', 'ADH1B', 'ADH1C']
##    #paralog = [paralog[0]]
##    oldest_paralog = 'ADH1A'
##    Dis = 'None'
##    Dir = False
##    gBGC = False
##    Force = {4:0.0}
##    print(oldest_paralog)
##
##    #for oldest_paralog in paralog:
##    test_tri = TriGeneconv( newicktree, alignment_file, save_path, oldest_paralog, Force = False, Dis = Dis, Dir = Dir, gBGC = gBGC)#, nnsites = nsites)
##    x_rates = np.log([0.00030354302771, # N0_N1
##                      0.0348850389125,   # N2_Ba
##                      0.0613673848999,   # N1_N2
##                      0.0131843207937,    # N3_Or
##                      0.010962330683,   # N2_N3
##                      0.00496736099364,  # N4_Go
##                      0.00784978378118,  # N3_N4
##                      0.000933915293812,  # N5_Bo
##                      0.00297723023713,  # N4_N5
##                      0.00656242294756,  # N6_Hu
##                      8.05308221277e-07,  # N5_N6
##                      0.00105781666618]) # N6_Ch
##    x = np.concatenate((x_process, x_rates))
##    test_tri.update_by_x(x)
##    print (test_tri._loglikelihood())
    #test_tri.get_mle(True, True)
#    test.get_individual_summary()
    


##    a = test.get_scene()
##    b = test_tri.get_scene()
##
##    process_a = a['process_definitions'][1]
##    process_b = b['process_definitions'][1]
##
##    process_a_dict = {str(list(process_a['row_states'][i]) + list(process_a['column_states'][i])):process_a['transition_rates'][i] for i in range(len(process_a['row_states']))}
##    process_b_dict = {str(list(process_b['row_states'][i]) + list(process_b['column_states'][i])):process_b['transition_rates'][i] for i in range(len(process_b['row_states']))}
##
##    print(process_a_dict == process_b_dict)
##
##    process_a = a['process_definitions'][0]
##    process_b = b['process_definitions'][2]
##
##    process_a_dict = {str(list(process_a['row_states'][i]) + list(process_a['column_states'][i])):process_a['transition_rates'][i] for i in range(len(process_a['row_states']))}
##    process_b_dict = {str(list(process_b['row_states'][i]) + list(process_b['column_states'][i])):process_b['transition_rates'][i] for i in range(len(process_b['row_states']))}
##
##    print(process_a_dict == process_b_dict)
##
##    tree_b_dict = {'_'.join([str(b['tree']['row_nodes'][i]), str(b['tree']['column_nodes'][i])]):b['tree']['edge_rate_scaling_factors'][i] for i in range(len(b['tree']['row_nodes']))}
##
##    for it in range(len(a['tree']['row_nodes'])):
##        row_node_name = test.tree.num_to_node[a['tree']['row_nodes'][it]]
##        col_node_name = test.tree.num_to_node[a['tree']['column_nodes'][it]]
##        if row_node_name[0] == 'D':
##            translated_row_name = 'N' + str(int(row_node_name[1:]) - 1)
##        elif row_node_name[0] == 'N':
##            translated_row_name = 'N' + str(int(row_node_name[1:]) + 2)
##        else:
##            print('check')
##
##        if col_node_name[0] == 'D':
##            translated_col_name = 'N' + str(int(col_node_name[1:]) - 1)
##        elif col_node_name[0] == 'N':
##            translated_col_name = 'N' + str(int(col_node_name[1:]) + 2)
##        else:
##            translated_col_name = col_node_name[:2] + '_'
##
##        tr_row_num = test_tri.node_to_num[translated_row_name]
##        tr_col_num = test_tri.node_to_num[translated_col_name]
##        tr_rate = tree_b_dict[str(tr_row_num) + '_' + str(tr_col_num)]
##
##        print(row_node_name, col_node_name, a['tree']['edge_rate_scaling_factors'][it],
##              tr_rate, a['tree']['edge_rate_scaling_factors'][it] == tr_rate,
##              test_tri.edge_to_blen[(translated_row_name, translated_col_name)])


#######################
    ###########Class 1 ADH Genes
#######################
##    gene_to_orlg_file = '../test/ADH1GeneToOrlg.txt'
##    #alignment_file = '../../ADH1Genes/Alignment/Intron_5_all/Intron_5_all_Mafft_gap_removed.fasta'
##    alignment_file = '../../ADH1Genes/Alignment/Concatenated_all_exon.fasta'
##    save_file = '../test/save/ADH1_all_exon_save.txt'
##    
##    data = Data(alignment_file, gene_to_orlg_file)
##    
##    tree_newick = '../test/PrimateTest.newick'
##    DupLosList = '../test/PrimateTestDupLost.txt'
##    tree = Tree(tree_newick, DupLosList)
##
##    node_to_pos = {'D1':0, 'D2':0, 'D3':1, 'D4':2, 'L1':2}
##    terminal_node_list = ['Chinese_Tree_Shrew', 'Macaque', 'Olive_Baboon', 'Orangutan', 'Gorilla', 'Human']
##    tree.get_configurations(terminal_node_list, node_to_pos)
##
##    pm_model = 'HKY'
##    x_js = np.log([0.5, 0.5, 0.5, 9.5, 4.9])
##    n_orlg = 9
##    IGC_pm = 'One rate'
##    n_js = 5
##    jsmodel = JSModel(n_js, x_js, pm_model, n_orlg, IGC_pm)
##    nsites = 5
##
##    test = JSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList, n_js, x_js, pm_model, n_orlg, IGC_pm,
##                      node_to_pos, terminal_node_list, save_file, nsites = nsites)
##    self = test
##    #process_definitions, conf_list = get_process_definitions(self.tree, self.jsmodel)
##    
##    conf_list = count_process(test.tree.node_to_conf)
##    configuration = conf_list[0]
###    process = jsmodel.get_process_definition(configuration)
##
##    edge_derivative = False
##    data = test.data
##    tree = test.tree
##    nsites = test.data.nsites
##    
##    data_type = test.jsmodel.PMModel.data_type
##    
##    observable_nodes, observable_axes, iid_observations = get_iid_observations(test.data, test.tree, test.data.nsites, test.jsmodel.PMModel.data_type)
##    
###    print test._loglikelihood()
##    test.get_mle(True, True)



