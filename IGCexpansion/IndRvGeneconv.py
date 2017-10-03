# A lazy compromise to implement HKY + rate variation model with IS-IGC model
# Xiang Ji
# xji3@ncsu.edu
from __future__ import print_function, absolute_import
from IndCodonGeneconv import IndCodonGeneconv
import numpy as np
import scipy
from copy import deepcopy
from functools import partial

class IndRvGeneconv:
    supported_model = ['HKY']
    
    def __init__(self, tree_newick, alignment_files, paralog, save_names, x,
                 Model = 'HKY',
                 nnsites = None, clock = False, Force = None,
                 post_dup = 'N1'):
        assert(len(alignment_files) == 3)
        assert(Model in IndRvGeneconv.supported_model)
        self.IndGeneconv_list = [IndCodonGeneconv( newicktree, alignment_files[alignment_iter], paralog, Model = model, Force = Force, clock = clock, save_name = save_names[alignment_iter]) \
                                 for alignment_iter in range(3)]
        self.model     = Model
        self.x         = np.array(x)
        # self.x = self.x_process + np.log([r2, r3]) + self.x_rates

        self.x_rates   = None
        self.x_process = None
        self.x_list    = None

        self.auto_save = 0

        self.initialize()

    def initialize(self):
        self.unpack_x(self.x)

    def update_by_x(self, x = None):
        if x is None:
            self.unpack_x(self.x)
        else:
            self.unpack_x(x)

    def unpack_x(self, x):
        assert(len(x) == len(self.x))
        self.x = x

        if self.model == 'HKY':
            self.x_process = np.concatenate((self.x[:4], [self.x[6]]))
            self.x_rates = self.x[7:]

            # Now construct x arrays
            self.x_list = [np.concatenate((self.x_process, self.x_rates + i)) for i in [0.0, self.x[4], self.x[5]]]
            for i in range(3):
                self.IndGeneconv_list[i].update_by_x(self.x_list[i])
                #self.IndGeneconv_list[i]._loglikelihood()
        else:
            raise Exception('Model has not been implemented!')

    def _loglikelihood(self, store = True, edge_derivative = False):
        self.unpack_x(self.x)
        ll =0.0
        if edge_derivative:
            div = np.array([0.0]*len(self.x_rates))
        else:
            div = []
        for i in range(3):
            site_ll, site_div = self.IndGeneconv_list[i]._loglikelihood(store, edge_derivative)
            ll += site_ll
            div += np.array(site_div)
        return ll, div

    def loglikelihood_and_gradient(self, display = False):
        self.update_by_x()
        delta = 1e-8
        x = deepcopy(self.x)  # store the current x array

        fn = self._loglikelihood

        ll, edge_derivs = fn(edge_derivative = True)
        
        m = len(self.x) - len(self.x_rates)

        # use finite differences to estimate derivatives with respect to these parameters
        other_derivs = []
        
        for i in range(m):
            x_plus_delta = np.array(self.x)
            x_plus_delta[i] += delta
            self.update_by_x(x_plus_delta)
            ll_delta, _ = fn(store = True, edge_derivative = False)
            d_estimate = (ll_delta - ll) / delta           
            other_derivs.append(d_estimate)
            # restore self.x
            self.update_by_x(deepcopy(x))
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
        if self.auto_save == 2:
            self.save_x()
            self.auto_save = 0
        return f, g

    def objective_wo_derivative(self, display, x):
        self.update_by_x(x)
        ll = self._loglikelihood()[0]

        if display:
            print ('log likelihood = ', ll)
            if self.clock:
                print ('Current x_clock array = ', self.x_clock)
            else:
                print ('Current x array = ', self.x)

        return -ll

    def get_mle(self, display = True, derivative = True):
        bnds = [(None, 0.0)] * 3 + [(None, None)] * (len(self.x) - 3)
        if derivative:
            f = partial(self.objective_and_gradient, display)
        else:
            f = partial(self.objective_wo_derivative, display)
        guess_x = self.x

        if derivative:
            result = scipy.optimize.minimize(f, guess_x, jac = True, method = 'L-BFGS-B', bounds = bnds)
        else:
            result = scipy.optimize.minimize(f, guess_x, jac = False, method = 'L-BFGS-B', bounds = bnds)

        print(result)
        return result

    def save_x(self):
        print
        
        
if __name__ == '__main__':
    
######################################################################################
######################################################################################
######################################################################################
    
    paralog = ['EDN', 'ECP']
    Force = None
    alignment_files = ['../test/EDN_ECP_Cleaned_CS_1.fasta', '../test/EDN_ECP_Cleaned_CS_2.fasta', '../test/EDN_ECP_Cleaned_CS_3.fasta']
    save_names = ['../test/EDN_ECP_HKY_HKY_rv_CS_' + str(i+1) + '_save.txt' for i in range(3)]
    newicktree = '../test/input_tree.newick'
    Force = None
    model = 'HKY'

    x = [-0.71492952, -0.55683374, -0.69445217,  0.74896689,  0.41813417, 0.44344036,  0.88867816,
         -2.56817753, -2.94618383, -3.26516553, -5.07067464, -4.81021362, -3.81058324, -5.68969739, -5.59586065]
    # This is the mle
    # log likelihood = -1713.0603332087344

    test = IndRvGeneconv( newicktree, alignment_files, paralog, save_names, x, Model = model, Force = Force, clock = None)
    self = test
    store = True
    display = True
    print (test._loglikelihood(True, True))
    print (test.loglikelihood_and_gradient(True))
    print (test.get_mle())
