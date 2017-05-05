# This script is used to infer IGC tract using HMM
# Xiang Ji
# xji3@ncsu.edu

import numpy as np
import scipy, scipy.optimize, scipy.linalg
from scipy.misc import logsumexp
from scipy.linalg import expm
from functools import partial
from math import floor
import os, sys

class HMMTract:
    def __init__(self, IGC_sitewise_lnL_file, Force_sitewise_lnL_file,
                 State_List, Total_blen, tau, seq_index_file):
        self.IGC_sitewise_lnL   = self.read_lnL(IGC_sitewise_lnL_file)
        self.Force_sitewise_lnL = self.read_lnL(Force_sitewise_lnL_file)
        self.StateList          = State_List
        self.L                  = Total_blen
        self.seq_index          = self.read_seq_index_file(seq_index_file)       

        # Now inferrence related parameters
        self.tau = tau         # estimated Tau value from MG94+IGC independent site model
        # Ptr is a function of the neighboring codon's distance on the chromosome
        # For example, when two codons have an intron between them, the Ptr is different
        # Decided not to have this stored
#        self.Ptr = None        # Transition probability matrix
        self.Emi = None        # Emission probability matrix
        self.eta = None        # IGC initiation rate
        self.tract_p = None    # IGC tract distribution p as in Geo(p)

        self.x      = None        # log array to store eta and tract_p values
        self.is_mle = False

        self.Forward_mat  = None   # To store values in Forward algorithm
        self.Backward_mat = None   # To store values in Backward algorithm

        self.init_parameters()

    def read_seq_index_file(self, seq_index_file):
        # The index should have columns:
        # nt_index, codon #, codon site for coding sequence
        seq_index = np.loadtxt(seq_index_file, dtype = int)
        # Dimension should match
        assert(seq_index.shape[0] == len(self.IGC_sitewise_lnL) * 3)
        return seq_index

    def init_parameters(self):
        self.update_by_x(np.log([self.tau * 0.05, 0.05]))
        assert(len(self.IGC_sitewise_lnL) == len(self.Force_sitewise_lnL))

    def read_lnL(self, sitewise_lnL_file):
        assert(os.path.isfile(sitewise_lnL_file))
        pos = []
        ll  = []
        with open(sitewise_lnL_file, 'rb') as f:
            for line in f:
                if line[0] == '#':
                    continue
                items = line.replace('\n', '').split('\t')
                pos.append(int(items[0]))
                ll.append(float(items[1]))

        return ll

    def update_by_x(self, x):
        assert(len(x) == 2)
        self.x = x
        self.eta = np.exp(x[0])
        self.tract_p = np.exp(x[1])
        self.get_Ptr_analytical()
        self.get_Emi()
        
    def get_marginal_state_distn(self):
        assert( self.eta != None and self.tract_p != None)
        P_S0 = np.exp( -2.0 * self.eta / self.tract_p * self.L)
        P_S1 = 1.0 - P_S0
        return [P_S0, P_S1]

    def get_Ptr(self):
        p = self.tract_p
        Q = self.eta * np.matrix([[-1.0 - 1.0/p, 1.0,   1.0, 1.0/p - 1.0],
                                  [0.0, -1.0/p, 0.0, 1.0/p],
                                  [0.0, 0.0, -1.0/p, 1.0/p],
                                  [0.0, 0.0, 0.0, 0.0]], dtype = float)
        P_mat = scipy.linalg.expm2(2 * Q * self.L)
        
        distn = self.get_marginal_state_distn()
        Ptr = np.matrix([[P_mat[0, 0] / distn[0], P_mat[0, 1] / distn[0] ],
                         [P_mat[0, 2] / distn[1], P_mat[0, 3] / distn[1] ]])

        #self.Ptr = np.log(Ptr)
        return np.log(Ptr)

    def get_Ptr_n(self, n):
        p = self.tract_p
        Q = self.eta * np.matrix([[-(2.0 - (1-p)**n)/p, (1.0 - (1-p)**n)/p, (1.0 - (1-p)**n)/p, (1 - p)**n/p],
                                  [0.0, -1.0/p, 0.0, 1.0/p],
                                  [0.0, 0.0, -1.0/p, 1.0/p],
                                  [0.0, 0.0, 0.0, 0.0]], dtype = float)
        P_mat = scipy.linalg.expm2(2 * Q * self.L)
        
        distn = self.get_marginal_state_distn()
        Ptr = np.matrix([[P_mat[0, 0] / distn[0], P_mat[0, 1] / distn[0] ],
                         [P_mat[0, 2] / distn[1], P_mat[0, 3] / distn[1] ]])

        return np.log(Ptr)

#        self.Ptr = np.log(Ptr)

    def get_Ptr_n_analytical(self, n):
        etl = np.exp(-2 * self.L * self.eta / self.tract_p)
        eel = np.exp(-2 * self.L * self.eta * (1.0 - (1.0 - self.tract_p)**n) / self.tract_p)
        Ptr = np.matrix([[eel, 1.0 - eel], [etl*(1.0 - eel)/(1.0 - etl), (1.0 - etl - (1.0 - eel) * etl)/(1.0 - etl)]], dtype = float)

        return np.log(Ptr)
#        self.Ptr = np.log(Ptr)

    def get_Ptr_analytical(self):
        etl = np.exp(-2 * self.L * self.eta / self.tract_p)
        eel = np.exp(-2 * self.L * self.eta)
        Ptr = np.matrix([[eel, 1.0 - eel], [etl*(1.0 - eel)/(1.0 - etl), (1.0 - etl - (1.0 - eel) * etl)/(1.0 - etl)]], dtype = float)
        #self.Ptr = np.log(Ptr)
        return np.log(Ptr)

    def get_Emi(self):
        self.Emi = np.zeros((len(self.StateList), len(self.IGC_sitewise_lnL)), dtype = float)
        distn = self.get_marginal_state_distn()
        for i in range(len(self.IGC_sitewise_lnL)):
            emission_0 = self.Force_sitewise_lnL[i]
            emission_1 = self.IGC_sitewise_lnL[i] + np.log(1.0 - np.exp(emission_0 - self.IGC_sitewise_lnL[i]) * distn[0]) - np.log(distn[1])

            if 1.0 - np.exp(emission_0 - self.IGC_sitewise_lnL[i]) * distn[0] < 0:
                print self.x
                sys.exit('something is wrong')

            self.Emi[:, i] = np.array([emission_0, emission_1])
            
        

    def Forward(self, display, x): # lnL by forward algorithm
        # update parameters first
        self.update_by_x(x)

        distn = self.get_marginal_state_distn()

        # Now create a 2 by nsites array for the dynamic programing
        lnL_array = np.zeros((len(self.StateList), len(self.IGC_sitewise_lnL)), dtype = float)

        # Now add in initial distribution
        lnL_array[:, 0] = np.log(distn) + self.Emi[:, 0]



        # Now do the forward step
        for i in range(len(self.IGC_sitewise_lnL) - 1):
            # Now calculate ln transition probabilities
            n = floor((self.seq_index[3*i + 3, 0] - self.seq_index[3*i, 0]) / 3)
            Ptr = self.get_Ptr_n_analytical(n)
            
            emission_0 = self.Emi[0, i + 1]
            emission_1 = self.Emi[1, i + 1]
            
            new_cond_lnL_0 = emission_0 + logsumexp([lnL_array[0, i] + Ptr[0, 0], lnL_array[1, i] + Ptr[1, 0]])
            new_cond_lnL_1 = emission_1 + logsumexp([lnL_array[0, i] + Ptr[0, 1], lnL_array[1, i] + Ptr[1, 1]])
            lnL_array[:, i + 1] = np.array([new_cond_lnL_0, new_cond_lnL_1])

        self.Forward_mat = lnL_array
        ll = lnL_array[0, -1] + np.log(sum(np.exp(lnL_array[:, -1] - lnL_array[0, -1])))
        if display:
            print '\t'.join([ str(item) for item in [ll] + list(np.exp(self.x))])

        return -ll

    def objective_1D(self, display, ln_p):
        x = np.array([np.log(self.tau) + ln_p[0], ln_p[0]])
        self.update_by_x(x)
        return self.Forward(display, x)
    

    def get_mle(self, display = True, OneDimension = True):
        self.update_by_x(self.x)
        if not OneDimension:
            f = partial(self.Forward, display)
            guess_x = self.x
            bnds = [(None, None), (None, 0.0)]

        else:
            f = partial(self.objective_1D, display)
            guess_x = self.x[1]
            bnds = [(None, 0.0)]

        result = scipy.optimize.minimize(f, guess_x, jac = False, method = 'L-BFGS-B', bounds = bnds)

        if display:
            print(result)
        if result['success']:
            self.is_mle = True
        return result

    def Viterbi(self):
        assert(self.is_mle)
        self.update_by_x(self.x)
        # same as the setup in Forward algorithm
        distn = self.get_marginal_state_distn()

        # Now create a 2 by nsites array for the Viterbi algorithm
        lnL_array = np.zeros((len(self.StateList), len(self.IGC_sitewise_lnL)), dtype = float)
        state_array = [[], []]

        # Now add in initial distribution
        lnL_array[:, 0] = np.log(distn) + self.Emi[:, 0]

        # Now do the Viterbi algorithm
        for i in range(len(self.IGC_sitewise_lnL) - 1):
            # Now calculate ln transition probabilities
            n = floor((self.seq_index[3*i + 3, 0] - self.seq_index[3*i, 0]) / 3)
            Ptr = self.get_Ptr_n_analytical(n)
            
            emission_0 = self.Emi[0, i + 1]
            emission_1 = self.Emi[1, i + 1]
            
            new_cond_lnL_0_list = [emission_0 + lnL_array[0, i] + Ptr[0, 0], emission_0 + lnL_array[1, i] + Ptr[1, 0]]
            lnL_array[0, i + 1] = max(new_cond_lnL_0_list)
            new_state = new_cond_lnL_0_list.index(max(new_cond_lnL_0_list))
            new_state_array_0 = state_array[new_state] + [new_state]

            new_cond_lnL_1_list = [emission_1 + lnL_array[0, i] + Ptr[0, 1], emission_1 + lnL_array[1, i] + Ptr[1, 1]]
            lnL_array[1, i + 1] = max(new_cond_lnL_1_list)
            new_state = new_cond_lnL_1_list.index(max(new_cond_lnL_1_list))
            new_state_array_1 = state_array[new_state] + [new_state]

            state_array = [new_state_array_0, new_state_array_1]

        last_state = list(lnL_array[:, -1]).index(max(lnL_array[:, -1]))
        Viterbi_path = state_array[last_state] + [last_state]

##        if display:
##            print

        return lnL_array, Viterbi_path
        

    def Backward(self):
        assert(self.is_mle)

        # S state distribution
        distn = self.get_marginal_state_distn()

        # create a 2 by nsites array for Backward algorithm storage
        # it will be passed to self.Backward_mat later
        lnL_array = np.zeros((len(self.StateList), len(self.IGC_sitewise_lnL)), dtype = float)

        # Now add in initial distribution
        n = floor((self.seq_index[-3, 0] - self.seq_index[-6, 0]) / 3)
        Ptr = self.get_Ptr_n_analytical(n)
        lnL_array[0, -1] = logsumexp([self.Emi[0, -1] + Ptr[0, 0], self.Emi[1, -1] + Ptr[0, 1]])
        lnL_array[1, -1] = logsumexp([self.Emi[0, -1] + Ptr[1, 0], self.Emi[1, -1] + Ptr[1, 1]])
        # Now do the backward steps
        for i in range(len(self.IGC_sitewise_lnL) - 1):
            # Now calculate ln transition probabilities
            n = floor((self.seq_index[-3*i - 3, 0] - self.seq_index[-3*i - 6, 0]) / 3)
            Ptr = self.get_Ptr_n_analytical(n)
            
            emission_0 = self.Emi[0, -(2 + i)]
            emission_1 = self.Emi[1, -(2 + i)]

            new_cond_lnL_0 = logsumexp([emission_0 + Ptr[0, 0] + lnL_array[0, -(i + 1)], \
                                        emission_1 + Ptr[0, 1] + lnL_array[1, -(i + 1)]])
            new_cond_lnL_1 = logsumexp([emission_0 + Ptr[1, 0] + lnL_array[0, -(i + 1)], \
                                        emission_1 + Ptr[1, 1] + lnL_array[1, -(i + 1)]])
            
            lnL_array[:, -(i + 2)] = np.array([new_cond_lnL_0, new_cond_lnL_1])
        self.Backward_mat = lnL_array
        return lnL_array

    def get_posterior(self):
        self.Backward()

        # initialize lnL array
        lnL_array = np.zeros((len(self.StateList), len(self.IGC_sitewise_lnL)), dtype = float)

        # Fill in the end position first
        lnL_array[:, -1] = self.Forward_mat[:, -1]

        # Fill in other positions
        lnL_array[:, :-1] = self.Backward_mat[:, 1:] + self.Forward_mat[:, :-1]
        lnL_array[0, :]   = lnL_array[0, :] - self.IGC_sitewise_lnL
        lnL_array[1, :]   = lnL_array[1, :] - self.IGC_sitewise_lnL

        return lnL_array
        
