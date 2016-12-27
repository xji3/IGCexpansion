# A separate file for JSLikelihood
# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Xiang Ji
# xji3@ncsu.edu

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

class JSGeneconv:
    def __init__(self, alignment_file, gene_to_orlg_file, # Data input
                 tree_newick, DupLosList,                 # Tree input
                 n_js, x_js, pm_model, n_orlg, IGC_pm,    # JSModel input
                 node_to_pos, terminal_node_list,         # Configuration input
                 save_file,                               # Auto save file
                 root_by_dup = False):   # JSModel input

        
        self.tree = Tree(tree_newick, DupLosList)
        self.data = Data(alignment_file, gene_to_orlg_file)
        self.jsmodel = JSModel(n_js, x_js, pm_model, n_orlg, IGC_pm)
        self.node_to_pos = node_to_pos
        self.terminal_node_list = terminal_node_list
        self.root_by_dup = root_by_dup
        self.save_file = save_file
        if os.path.isfile(self.save_file):
            self.initialize_by_save()
            print ('Loaded paramters from ' + self.save_file)
        else:            
            if self.root_by_dup:
                self.x = np.concatenate((np.log([0.1] * len(self.tree.edge_list)), x_js))
            else:
                self.x = np.concatenate((np.log([0.1] * (len(self.tree.edge_list) - 1)), x_js))

        self.tree.get_configurations(self.terminal_node_list, self.node_to_pos)
        assert(self.jsmodel.n_orlg == self.tree.n_orlg)  # make sure n_orlg got updated
        self.ll   = None    # used to store log likelihood
        self.auto_save = 0  # used to control auto save frequency
        
    def unpack_x(self, x):
        self.x = x
        if self.root_by_dup:
            x_rate = x[:len(self.tree.edge_list)]
            x_js   = x[len(self.tree.edge_list):]
        else:
            x_rate = x[:(len(self.tree.edge_list) - 1)]
            x_js   = x[(len(self.tree.edge_list) - 1):]
        self.unpack_x_rates(x_rate)
        self.jsmodel.unpack_x_js(x_js)

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

    def get_scene(self, nsites = None):
        state_space_shape = self.jsmodel.state_space_shape
        process_definitions, conf_list = get_process_definitions(self.tree, self.jsmodel)
        self.tree.get_tree_process(conf_list)
        if not self.root_by_dup:
            prior_feasible_states, distn = self.jsmodel.get_prior(self.tree.node_to_conf['N0'])
        else:
            print 'This case is not implemented yet. Check get_scene function in JSGeneconv class.'

        if nsites is None:
            observable_nodes, observable_axes, iid_observations = get_iid_observations(self.data, self.tree, self.data.nsites, self.jsmodel.PMModel.data_type)
        else:
            observable_nodes, observable_axes, iid_observations = get_iid_observations(self.data, self.tree, nsites, self.jsmodel.PMModel.data_type)
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

    def _loglikelihood(self, nsites = None, edge_derivative = False):
        if nsites is None:
            scene = self.get_scene()
        else:
            scene = self.get_scene(nsites)
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

    
    def loglikelihood_and_gradient(self, nsites = None, display = False):
        delta = 1e-8
        x = deepcopy(self.x)
        ll, edge_derivs = self._loglikelihood(nsites = nsites, edge_derivative = True)
        if self.root_by_dup:
            x_rate_derivs = {self.tree.edge_list[i]:edge_derivs[i] for i in range(len(edge_derivs))}
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
            x_plus_delta = np.array(self.x)
            x_plus_delta[i] += delta
            self.unpack_x(x_plus_delta)
            ll_delta, _ = self._loglikelihood(nsites, False)
            d_estimate = (ll_delta - ll) / delta
            other_derivs.append(d_estimate)
            # restore self.x
            self.unpack_x(x)
        other_derivs = np.array(other_derivs)
        self.ll = ll
        f = -ll
        g = -np.concatenate((x_rate_derivs, other_derivs))
        if display:
            print ('log likelihood = ', ll)
            print ('Edge derivatives = ', x_rate_derivs)
            print ('other derivatives:', other_derivs)
            print ('Current x array = ', self.x)
            
        return f, g

    def objective_and_gradient(self,display, x):
        self.unpack_x(x)
        f, g = self.loglikelihood_and_gradient(display = display)
        self.auto_save += 1
        if self.auto_save == 5:
            self.save_x()
            self.auto_save = 0
        return f, g

    def objective_wo_gradient(self, display, x):
        self.unpack_x(x)
        ll = self._loglikelihood()[0]
        if display:
            print ('log likelihood = ', ll)
            print ('Current x array = ', self.x)

        return -ll
            
    def get_mle(self, display = True, derivative = True):
        self.unpack_x(self.x)  # do one more update first
        if derivative:
            f = partial(self.objective_and_gradient, display)
        else:
            f = partial(self.objective_wo_derivative, display)

        guess_x = self.x
        if self.root_by_dup:
            bnds  = [(None, None)] * len(self.tree.edge_list)
        else:
            bnds  = [(None, None)] * (len(self.tree.edge_list) - 1)

        bnds.extend([(None, -0.001)] * 3)
        bnds.extend([(None, None)] * (len(self.jsmodel.x_js) - 3))

        if derivative:
            result = scipy.optimize.minimize(f, guess_x, jac = True, method = 'L-BFGS-B', bounds = bnds)
        else:
            result = scipy.optimize.minimize(f, guess_x, jac = False, method = 'L-BFGS-B', bounds = bnds)

        print(result)
        return result

    def save_x(self):
        save = self.x
        np.savetxt(open(self.save_file, 'w+'), save.T)

    def initialize_by_save(self):
        self.x = np.loadtxt(open(self.save_file, 'r'))
        self.unpack_x(self.x)
        
                    
        


if __name__ == '__main__':

#######################
    ###########Yeast
#######################
##    gene_to_orlg_file = '../test/YLR406_YDL075W_GeneToOrlg.txt'
##    alignment_file = '../test/YLR406C_YDL075W_alignment.fasta'
##
##    tree_newick = '../test/YeastTree.newick'
##    DupLosList = '../test/YeastTestDupLost.txt'
##    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
##    node_to_pos = {'D1':0}
##
##    pm_model = 'HKY'
##    x_js = np.log([ 0.49978115,   0.60208772,   0.41240341,  10.35588244,   8.01054376])
##    n_orlg = 3
##    IGC_pm = 'One rate'
##    n_js = 2
##
##    test = JSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList, n_js, x_js, pm_model, n_orlg, IGC_pm,
##                      node_to_pos, terminal_node_list)
##
##    self = test
###    print test.get_scene(100)
##    print test._loglikelihood(edge_derivative = True)
##    print test.loglikelihood_and_gradient()
##    test.get_mle()
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
##    tree = Tree(tree_newick, DupLosList)
##
##    node_to_pos = {'D1':0, 'D2':0}
##    terminal_node_list = ['Baboon', 'Orangutan', 'Gorilla', 'Bonobo', 'Chimpanzee', 'Human']
##    tree.get_configurations(terminal_node_list, node_to_pos)
##
##    pm_model = 'HKY'
##    x_js = np.log([ 0.49978115,   0.60208772,   0.41240341,  10.35588244,   8.01054376])
##    n_orlg = 5
##    IGC_pm = 'One rate'
##    n_js = 3
##    jsmodel = JSModel(n_js, x_js, pm_model, n_orlg, IGC_pm)
##
##    test = JSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList, n_js, x_js, pm_model, n_orlg, IGC_pm,
##                      node_to_pos, terminal_node_list)
##    self = test
##    x = np.concatenate((np.log([0.2] * (len(test.tree.edge_list) - 1)), x_js))
##    test.unpack_x(x)
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
##    data_type = test.jsmodel.PMModel.data_type
##    
##    observable_nodes, observable_axes, iid_observations = get_iid_observations(test.data, test.tree, test.data.nsites, test.jsmodel.PMModel.data_type)
##    print test._loglikelihood()
##    #test.get_mle()


#######################
    ###########Class 1 ADH Genes
#######################
    gene_to_orlg_file = '../test/ADH1GeneToOrlg.txt'
    alignment_file = '../../ADH1Genes/Alignment/Intron_5_all/Intron_5_all_Mafft_gap_removed.fasta'
    save_file = '../test/save/ADH1_intron_5_save.txt'
    
    data = Data(alignment_file, gene_to_orlg_file)
    
    tree_newick = '../test/PrimateTest.newick'
    DupLosList = '../test/PrimateTestDupLost.txt'
    tree = Tree(tree_newick, DupLosList)

    node_to_pos = {'D1':0, 'D2':0, 'D3':1, 'D4':2, 'L1':2}
    terminal_node_list = ['Chinese_Tree_Shrew', 'Macaque', 'Olive_Baboon', 'Orangutan', 'Gorilla', 'Human']
    tree.get_configurations(terminal_node_list, node_to_pos)

    pm_model = 'HKY'
    x_js = np.log([0.3, 0.5, 0.2, 9.5, 4.9])
    n_orlg = 9
    IGC_pm = 'One rate'
    n_js = 5
    jsmodel = JSModel(n_js, x_js, pm_model, n_orlg, IGC_pm)

    test = JSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList, n_js, x_js, pm_model, n_orlg, IGC_pm,
                      node_to_pos, terminal_node_list, save_file)
    self = test
    #process_definitions, conf_list = get_process_definitions(self.tree, self.jsmodel)
    
    conf_list = count_process(test.tree.node_to_conf)
    configuration = conf_list[0]
#    process = jsmodel.get_process_definition(configuration)

    edge_derivative = False
    data = test.data
    tree = test.tree
    nsites = test.data.nsites
    data_type = test.jsmodel.PMModel.data_type
    
    observable_nodes, observable_axes, iid_observations = get_iid_observations(test.data, test.tree, test.data.nsites, test.jsmodel.PMModel.data_type)
    
#    print test._loglikelihood()
#    test.get_mle()



