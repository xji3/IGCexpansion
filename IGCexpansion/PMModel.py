# A separate class to represent Point mutation model
# Xiang Ji
# xji3@ncsu.edu
import numpy as np
import sys
from copy import deepcopy
from operator import mul
from itertools import product

bases = 'tcag'.upper()
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = [amino_acid for amino_acid in 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']

class PMModel:
    supported = ['HKY', 'MG94']                   # supported models
    # Constants for Sequence operations
    codon_table = dict(zip(codons, amino_acids))
    codon_nonstop = [codons[i] for i in range(len(codons)) if amino_acids[i] != '*']
    
    def __init__(self, model_name, x_pm, rate_variation, force = None):
        self.name           = model_name          # name of supported models
        self.x_pm           = x_pm                # an array of log() values
        self.data_type      = None                # used for get_iid_observation function in Func.py
        self.force          = force               # Used for parameter value constraint

        self.Q_mut          = None                # Point mutation Q matrix
        self.parameters     = dict()              # a dictionary to store all model parameters
        self.parameter_list = None
        self.rate_variation = rate_variation      # a bool indicator of whether rate variation among codon sites is considered

        self.normalizing_factor = None            # normalizing factor used in stationary distribution, should be 1 for HKY
        self.expected_rate      = None            # normalizing factor for rate matrix Q unit
        self.init_Q()


    def init_Q(self):
        assert(self.name in self.supported) # check if PM model is implemented
        if self.name == 'HKY':
            if self.rate_variation:
                self.parameter_list = ['Pi_A', 'Pi_C', 'Pi_G', 'Pi_T', 'kappa', 'r2', 'r3']
                # r1, r2, r3 are three rates for three nt positions in a codon.
                # r1 is always 1.0 and is used for the branch length unit.
                # The name of the parameters are suggested by Jeff.
            else:
                self.parameter_list = ['Pi_A', 'Pi_C', 'Pi_G', 'Pi_T', 'kappa']

            self.init_HKY_Q()
            self.data_type = 'nt'
            
        elif self.name == 'MG94':
            assert(not self.rate_variation)
            self.parameter_list = ['Pi_A', 'Pi_C', 'Pi_G', 'Pi_T', 'kappa', 'omega']
            self.data_type = 'cd'
            self.init_MG94_Q()
            

    def init_HKY_Q(self):
        # This function initialize Q matrix for HKY85 point mutation model
        # first, check x_pm size
        # x_pm = np.log([%AG, %A, %C, kappa])
        # HKY model has 3 parameters for nt frequency, and 1 parameter kappa for transition/transversion ratio
        assert(self.check_x_pm())
        assert(self.check_force())

        self.unpack_frequency()
        if self.force == None or not 3 in self.force:
            kappa = np.exp(self.x_pm[3])
        else:
            kappa = self.force[3]
        self.parameters['kappa'] = kappa

        # Now add in rate variation
        if self.rate_variation:
            if self.force == None or not 4 in self.force:
                r2 = np.exp(self.x_pm[4])
            else:
                r2 = self.force[4]

            if self.force == None or not 5 in self.force:
                r3 = np.exp(self.x_pm[5])
            else:
                r3 = self.force[5]

            self.parameters['r2'] = r2
            self.parameters['r3'] = r3

        # In order of
        # ACGT   A=0, C=1, G=2, T=3
        Qbasic = np.array([
            [0, 1.0, kappa, 1.0],
            [1.0, 0, 1.0, kappa],
            [kappa, 1.0, 0, 1.0],
            [1.0, kappa, 1.0, 0],
            ]) * np.array([self.parameters['Pi_' + nt] for nt in 'ACGT'])
        # assume PM at stationary
        stationary_distn = [ self.parameters['Pi_' + nt] for nt in 'ACGT' ]
        assert(np.isclose(sum(stationary_distn), 1.0))
        self.normalizing_factor = sum(stationary_distn)
        
        expected_rate = np.dot(stationary_distn, Qbasic.sum(axis = 1))
        self.Q_mut = Qbasic / expected_rate
        self.expected_rate = expected_rate

    def init_MG94_Q(self):
        # This function initialize Q matrix for MG94 codon substitution model
        # first , check x_pm size
        # x_pm = np.log([%AG, %A, %C, kappa, omega])
        # MG94 model has 3 parameters for codon frequency, 1 parameter kappa for transition/transversion ratio and
        # 1 parameter omega for natural nonsynonimous/synonious ratio selection
        assert(self.check_x_pm())
        assert(self.check_force())

        self.unpack_frequency()
        # kappa
        if self.force == None or not 3 in self.force:
            kappa = np.exp(self.x_pm[3])
        else:
            kappa = self.force[3]
        self.parameters['kappa'] = kappa

        # omega
        if self.force == None or not 4 in self.force:
            omega = np.exp(self.x_pm[4])
        else:
            omega = self.force[4]
        self.parameters['omega'] = omega

        # In order of
        # PMModel.codon_nonstop
        Qbasic = np.zeros((len(PMModel.codon_nonstop), len(PMModel.codon_nonstop)), dtype = float)
        # Now assign off diagonal entries by looping
        for row_iter in range(len(PMModel.codon_nonstop)):
            state_from = PMModel.codon_nonstop[row_iter]
            for col_iter in range(len(PMModel.codon_nonstop)):
                if row_iter == col_iter:
                    continue
                state_to = PMModel.codon_nonstop[col_iter]
                diff_iter = [i for i in range(3) if state_from[i] != state_to[i]]
                if len(diff_iter) > 1:
                    continue
                else:
                    nt_to = state_to[diff_iter[0]]
                    Qbasic[row_iter, col_iter] = self.parameters['Pi_' + nt_to]
                    if not self.is_transversion(state_from, state_to):
                        Qbasic[row_iter, col_iter] = Qbasic[row_iter, col_iter] * self.parameters['kappa']
                    if self.is_nonsynonymous(state_from, state_to):
                        Qbasic[row_iter, col_iter] = Qbasic[row_iter, col_iter] * self.parameters['omega']

        # assume stationary
        stationary_distn = [ reduce(mul, [self.parameters['Pi_' + b]  for b in codon], 1) for codon in PMModel.codon_nonstop ]
        self.normalizing_factor = sum(stationary_distn)
        stationary_distn = np.array(stationary_distn) / self.normalizing_factor
        expected_rate = np.dot(stationary_distn, Qbasic.sum(axis = 1))
        self.Q_mut = Qbasic / expected_rate
        self.expected_rate = expected_rate
        
        
    def is_transversion(self, state_from, state_to):
        assert(state_from != state_to)
        if self.data_type == 'nt':
            assert(len(state_from) == len(state_to) == 1)
            nt_from = state_from
            nt_to = state_to
        elif self.data_type == 'cd':
            assert(len(state_from) == len(state_to) == 3)
            diff_iter = [i for i in range(3) if state_from[i] != state_to[i]]
            assert(len(diff_iter) == 1)
            nt_from = state_from[diff_iter[0]]
            nt_to   = state_to[diff_iter[0]]
        else:
            sys.exit('Data type not supported!')

        return not (set([nt_from, nt_to]) == set(['A', 'G']) or set([nt_from, nt_to]) == set(['C', 'T']))

    def is_nonsynonymous(self, state_from, state_to):
        assert(state_from != state_to)
        assert(self.data_type == 'cd')
        assert(len(state_from) == len(state_to) == 3)
        return PMModel.codon_table[state_from] != PMModel.codon_table[state_to]
        
        

    def check_x_pm(self):
        if self.name == 'HKY':
            if self.rate_variation:
                return len(self.x_pm) == 6
            else:
                return len(self.x_pm) == 4
        elif self.name == 'MG94':
            return len(self.x_pm) == 5
        else:
            sys.exit('The point mutation model has not been implemented!')
            return False # Just in case

    def check_force(self):
        if self.force == None:
            return True
        else:
            if self.name == 'HKY':
                if self.rate_variation:
                    return all([key < 6 for key in self.force])
                else:
                    return all([key < 4 for key in self.force])
            elif self.name == 'MG94':
                return all([key < 5 for key in self.force])
            else:
                sys.exit('The point mutation model has not been implemented!')
                return False # Just in case sys.exit is missed, so that assert still captures.


    def unpack_frequency(self):
    	#%AG, %A, %C
        x_process = np.exp(self.x_pm[:3])
        if not self.force == None:
            for i in range(3):
                if i in self.force:
                    x_process[i] = self.force[i]
        pi_a = x_process[0] * x_process[1]
        pi_c = (1 - x_process[0]) * x_process[2]
        pi_g = x_process[0] * (1 - x_process[1])
        pi_t = (1 - x_process[0]) * (1 - x_process[2])
        self.parameters['Pi_A'] = pi_a
        self.parameters['Pi_C'] = pi_c
        self.parameters['Pi_G'] = pi_g
        self.parameters['Pi_T'] = pi_t


    def update_by_x_pm(self, new_x_pm):
        assert(len(self.x_pm) == len(new_x_pm))
        self.x_pm = deepcopy(new_x_pm)
        self.init_Q()

    def get_stationary_distn(self, state):
        if self.name == 'HKY':
            return self.get_HKY_stationary_distn(state)
        elif self.name == 'MG94':
            return self.get_MG94_stationary_distn(state)

    def get_MG94_stationary_distn(self, state):
        assert( -1 < state < len(PMModel.codon_nonstop))
        codon = PMModel.codon_nonstop[state]
        stationary_rate = reduce(mul, [self.parameters['Pi_' + b]  for b in codon], 1) / self.normalizing_factor
        return stationary_rate
        
        
    def get_HKY_stationary_distn(self, state):
        assert(-1 < state < 4)
        # 0:A, 1:C, 2:G, 3:T
        return self.parameters['Pi_' + 'ACGT'[state]]

    def get_HKY_transition_rate(self, transition, codon_site):
        state_from, state_to = transition
        if codon_site == 1:
            return self.Q_mut[state_from, state_to]
        elif codon_site == 2:
            assert(self.rate_variation)
            return self.Q_mut[state_from, state_to] * self.parameters['r2']
        elif codon_site == 3:
            assert(self.rate_variation)
            return self.Q_mut[state_from, state_to] * self.parameters['r3']
        else:
            sys.exit('Codon_site out of 3 is not allowed!')

    def __str__(self): # overide for print function
        return 'Point mutation model: ' + self.name + '\n' + \
               'Rate variation: ' + str(self.rate_variation) + '\n' + \
               'Point mutation parameters: ' + ' '.join([item + ' '+ str(self.parameters[item]) for item in self.parameter_list]) + '\n'

if __name__ == '__main__':
##    test = PMModel('HKY', np.log([0.3, 0.5, 0.2, 9.5]), False, {2:0.0})
##    self = test
##    print test.Q_mut
##    test.update_by_x_pm(np.log([0.1, 0.9, 0.3, 11.0]))
##    print test.Q_mut
##
##    test = PMModel('HKY', np.log([0.3, 0.5, 0.2, 9.5, 0.4, 1.4]), True)
##    self = test
##    print test.Q_mut
##    test.update_by_x_pm(np.log([0.1, 0.9, 0.3, 11.0, 0.4, 1.5]))
##    print test.Q_mut

    test = PMModel('MG94', np.log([0.3, 0.5, 0.2, 9.5, 0.3]), False)
    self = test
    test.init_Q()
    
