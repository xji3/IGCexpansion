# A separate class to represent Pair Site Joint State IGC model (PS JS IGC models)
# PS JS IGC model = IGC model + Point mutation model
# Xiang Ji
# xji3@ncsu.edu
import sys
from PMModel import PMModel
from PSIGCModel import PSIGCModel
import numpy as np
import itertools
from copy import deepcopy
from operator import mul
from scipy.sparse import lil_matrix
import scipy.sparse.linalg
from Common import *

class PSJSModel:
    supported = ['One rate']
    # Consider only two paralog case
    def __init__(self, x_js, pm_model, n_orlg, IGC_pm, rate_variation = False, force = None, n_js = 2):
        self.n_js   = n_js            # number of contemporaneous paralog states considered on each branch
        self.x_js   = x_js            # a concatenated vector storing x_pm + x_IGC        
        self.x_pm   = None            # x_pm vector for PMModel
        self.x_IGC  = None            # x_IGC vector for IGCModel
        self.IGC_force = None         # parameter value constraints on IGC
        self.force  = force           # parameter value constraint
        self.IGC_pm = IGC_pm
        self.rate_variation = rate_variation # bool indicator of rate variation for point mutation model

        self.pm_model = pm_model      # name of point mutation model
        self.n_orlg   = n_orlg        # total number of ortholog groups

        self.PMModel  = None          # PMModel class instance for point mutation model
        # TODO: use PSIGCModel class
#        self.PSIGCModel = None        # PSIGCModel class instance for pair state IGC model

        self.state_space_shape = None # initialized in init_models() function

        self.init_models()


    def unpack_x_js(self, x_js):
        # first, check if the models are supported
        assert(self.pm_model in PMModel.supported)
        assert(self.IGC_pm in PSJSModel.supported)
        if self.pm_model == 'HKY':
            num_x_pm = 4
            if self.rate_variation:
                num_x_pm += 2
        else:
            sys.exit( 'The point mutation model is not supported.')

        if self.IGC_pm == 'One rate':
            num_x_IGC = 2
        else:
            sys.exit( 'The IGC parameterization has not been implemented.')
            
        self.x_pm = x_js[:num_x_pm]
        self.x_IGC = x_js[num_x_pm:]
        self.x_js = x_js
        assert(num_x_pm + num_x_IGC == len(self.x_js))

    def update_by_x_js(self, new_x_js):
        self.unpack_x_js(new_x_js)
        self.PMModel.update_by_x_pm(self.x_pm)
        #self.PSIGCModel.update_by_x_IGC(self.x_IGC)


    def divide_force(self):
        if self.force == None:
            return None, None
        # first, check if the models are supported
        assert(self.pm_model in PMModel.supported)
        assert(self.IGC_pm in PSJSModel.supported)
        if self.pm_model == 'HKY':
            num_x_pm = 4
            if self.rate_variation:
                num_x_pm += 2
        else:
            sys.exit('The point mutation model is not supported.')

        if self.IGC_pm == 'One rate':
            num_x_IGC = 2
        else:
            sys.exit( 'The IGC parameterization has not been implemented.')

        pm_force = dict()
        IGC_force = dict()
        for key in self.force:
            if key < num_x_pm:
                pm_force[key] = self.force[key]
            else:
                IGC_force[key - num_x_pm] = self.force[key]

        if not pm_force.keys():
            pm_force = None
        if not IGC_force.keys():
            IGC_force = None

        return pm_force, IGC_force

    def init_models(self):
        self.unpack_x_js(self.x_js)
        if self.pm_model == 'HKY':
            self.state_space_shape = [4,4] * self.n_js
        else:
            sys.exit('The point mutation model has not been implemented.')

        pm_force, IGC_force = self.divide_force()
        self.PMModel = PMModel(self.pm_model, self.x_pm, self.rate_variation, pm_force)
        #self.PSIGCModel = PSIGCModel(self.x_IGC, self.n_orlg, self.IGC_pm, force = IGC_force)
        self.IGC_force = IGC_force
        assert( len(set(self.state_space_shape)) == 1) # now consider only same state space model
        self.update_by_x_js(self.x_js)

    def is_transition_compatible(self, transition):
        # only consider two paralogs for now
        assert(len(self.state_space_shape) == 2 * self.n_js)
        assert(len(transition) == 2)
        state_from, state_to = transition
        assert(len(state_from) == len(state_to) == len(self.state_space_shape))
        if state_from == state_to:
            return False

        # state = (ia, ib, ja, jb) with two paralogs i, j and two positions a, b
        # Now get positions in state that are different in state_from, state_to
        pos_list = [i for i in range(len(state_from)) if state_from[i] != state_to[i]]
        if len(pos_list) > 2:
            return False
        elif len(pos_list) == 2: # only IGC can copy two sites over 
            if pos_list == [0, 1] and state_to[0] == state_from[2] and state_to[1] == state_from[3] \
               or pos_list == [2, 3] and state_to[2] == state_from[0] and state_to[3] == state_from[1]:
                return True
            else:
                return False
        elif len(pos_list) == 1: # one state change can be from mutation / IGC
            return True
        else:
            sys.exit('Check is_transition_compatible function in PSJSModel class.')
        
        
    def cal_IGC_transition_rate(self, transition, n, codon_site_pair, proportion = False):
        assert(all([ 0 < i < 4 for i in codon_site_pair]))
        # n is the distance between two sites
        # n = 1, 2, 3, ...
        # Transition should be compatible
        assert(self.is_transition_compatible(transition))
        
        # Now get the two states
        state_from, state_to = transition

        # Now extend the codon_site_pair to have compatible size as states
        extended_cdn_site_pair = codon_site_pair * 2
        assert(len(extended_cdn_site_pair) == len(state_from))

        # Get positions in the state that differ
        pos_list = [i for i in range(len(state_from)) if state_from[i] != state_to[i]]

        if self.IGC_pm == 'One rate':
            if self.IGC_force is None:
                IGC_init, IGC_p = np.exp(self.x_IGC)
            else:
                exp_x_IGC = np.exp(self.x_IGC)
                for key in self.IGC_force:
                    exp_x_IGC[key] = self.IGC_force[key]
                IGC_init, IGC_p = exp_x_IGC
            # Now calculate IGC rate
            IGC_0_not_n = IGC_init / IGC_p * (1 - (1 - IGC_p) ** n)
            IGC_0_and_n = IGC_init / IGC_p * (1 - IGC_p) ** n
            if len(pos_list) == 1:
                pos = pos_list[0]
                codon_site = extended_cdn_site_pair[pos]
                same_paralog_other_pos, other_paralog_same_pos, other_paralog_other_pos = self.get_other_pos(pos)
                q_ij = self.PMModel.get_HKY_transition_rate((state_from[pos], state_to[pos]), codon_site)
                if state_to[pos] == state_from[other_paralog_same_pos] \
                   and state_from[same_paralog_other_pos] == state_from[other_paralog_other_pos]:
                    q_IGC = IGC_0_not_n + IGC_0_and_n
                elif state_to[pos] == state_from[other_paralog_same_pos] \
                   and state_from[same_paralog_other_pos] != state_from[other_paralog_other_pos]:
                    q_IGC = IGC_0_not_n
                else:
                    q_IGC = 0.0
            elif len(pos_list) == 2:
                q_ij = 0.0
                if pos_list == [0, 1] and state_to[0] == state_from[2] and state_to[1] == state_from[3] \
               or pos_list == [2, 3] and state_to[2] == state_from[0] and state_to[3] == state_from[1]:
                    q_IGC = IGC_0_and_n
                else:
                    sys.exit('Transition not compatible!')
            else:
                sys.exit('Transition not compatible! 2')

            if proportion:
                return q_IGC / (q_ij + q_IGC)
            else:
                return q_ij + q_IGC

        else:
            sys.exit('Cal_IGC_transition_rate not implemented yet.')
    
    def get_IGC_transition_rates_BF(self, n, codon_site_pair, proportion = False):
        assert(all([ 0 < i < 4 for i in codon_site_pair]))
        # This function is only for checking code, it should not be used in real calculation
        for state_from in itertools.product(range(4), repeat = self.n_js * 2):
            for state_to in itertools.product(range(4), repeat = self.n_js * 2):
                if self.is_transition_compatible((state_from, state_to)):
                    yield state_from, state_to, self.cal_IGC_transition_rate((state_from, state_to), n, codon_site_pair, proportion)

    def get_IGC_transition_rates(self, n, codon_site_pair, proportion = False):
        assert(all([ 0 < i < 4 for i in codon_site_pair]))
        # This function is only for checking code, it should not be used in real calculation
        for state_from in itertools.product(range(4), repeat = self.n_js * 2):
            # Now visit only possible state_to
            # One single change
            for i in range(len(state_from)):
                for nt in range(4):
                    if nt == state_from[i]:
                        continue
                    new_state = list(deepcopy(state_from))
                    new_state[i] = nt
                    new_state = tuple(new_state)
                    yield state_from, new_state, self.cal_IGC_transition_rate((state_from, new_state), n, codon_site_pair, proportion)
            
            # two site changes
            if state_from[0] != state_from[2] and state_from[1] != state_from[3]:
                state_to = list(deepcopy(state_from))
                state_to[0] = state_from[2]
                state_to[1] = state_from[3]
                state_to = tuple(state_to)
                yield state_from, state_to, self.cal_IGC_transition_rate((state_from, state_to), n, codon_site_pair, proportion)

                state_to = list(deepcopy(state_from))
                state_to[2] = state_from[0]
                state_to[3] = state_from[1]
                state_to = tuple(state_to)
                yield state_from, state_to, self.cal_IGC_transition_rate((state_from, state_to), n, codon_site_pair, proportion)


    def is_compatible_pm_transition(self, transition):
        # This PM transition is for pre-duplication rate matrix construction
        # only consider two paralogs for now
        assert(len(self.state_space_shape) == 2*self.n_js)
        assert(len(transition) == 2)
        state_from, state_to = transition
        assert(len(state_from) == len(state_to) == len(self.state_space_shape))
        if state_from == state_to:
            return False

        # state = (ia, ib, ja, jb) with two paralogs i, j and two positions a, b
        # Now get positions in state that are different in state_from, state_to
        if state_from[0] != state_from[2] or state_from[1] != state_from[3] \
           or state_to[0] != state_to[2] or state_to[1] != state_to[3]:
            return False
        
        pos_list = [i for i in range(len(state_from) / 2) if state_from[i] != state_to[i]]
        if len(pos_list) == 1:
            return True
        else:
            return False

    def cal_PM_transition_rate(self, transition, codon_site_pair):
        assert(all([ 0 < i < 4 for i in codon_site_pair]))
        assert(self.is_compatible_pm_transition(transition))
        state_from, state_to = transition
        pos_list = [i for i in range(len(state_from) / 2) if state_from[i] != state_to[i]]
        pos = pos_list[0]
        codon_site = codon_site_pair[pos]
        return self.PMModel.get_HKY_transition_rate((state_from[pos], state_to[pos]), codon_site)
            
    def get_PM_transition_rates_BF(self, codon_site_pair):
        assert(all([ 0 < i < 4 for i in codon_site_pair]))
        # This function is only for checking code, it should not be used in real calculation
        for state_from in itertools.product(range(4), repeat = self.n_js):
            state_from = state_from + state_from
            for state_to in itertools.product(range(4), repeat = self.n_js):
                state_to = state_to + state_to
                if self.is_compatible_pm_transition((state_from, state_to)):
                    yield state_from, state_to, self.cal_PM_transition_rate((state_from, state_to), codon_site_pair) 

    def get_PM_transition_rates(self, codon_site_pair):
        assert(all([ 0 < i < 4 for i in codon_site_pair]))
        for state_from in itertools.product(range(4), repeat = self.n_js):
            state_from = state_from + state_from
            for i in range(2):
                for nt in range(4):
                    if state_from[i] == nt:
                        continue
                    state_to = list(deepcopy(state_from))
                    state_to[i] = state_to[i + 2] = nt
                    state_to = tuple(state_to)
                    yield state_from, state_to, self.cal_PM_transition_rate((state_from, state_to), codon_site_pair)
                     
        
    def get_IGC_process_definition(self, n, codon_site_pair = (1, 1), proportion = False):
        assert(all([ 0 < i < 4 for i in codon_site_pair]))
        row_states = []
        column_states = []
        transition_rates = []
        for row_state, col_state, transition_rate in self.get_IGC_transition_rates(n, codon_site_pair, proportion):
            row_states.append(row_state)#translate_four_nt_to_two_state(row_state))
            column_states.append(col_state)#translate_four_nt_to_two_state(col_state))
            transition_rates.append(transition_rate)

        if proportion:
            process_definition = dict(
                row_states = row_states,
                column_states = column_states,
                weights = transition_rates)
        else:
            process_definition = dict(
                row_states = row_states,
                column_states = column_states,
                transition_rates = transition_rates)
        return process_definition

    def get_PM_process_definition(self, codon_site_pair = (1, 1)):
        assert(all([ 0 < i < 4 for i in codon_site_pair]))
        row_states = []
        column_states = []
        transition_rates = []
        for row_state, col_state, transition_rate in self.get_PM_transition_rates(codon_site_pair):
            row_states.append(row_state)#translate_four_nt_to_two_state(row_state))
            column_states.append(col_state)#translate_four_nt_to_two_state(col_state))
            transition_rates.append(transition_rate)

        process_definition = dict(
            row_states = row_states,
            column_states = column_states,
            transition_rates = transition_rates)
        return process_definition
    
    def get_other_pos(self, pos):
        # return same_paralog_other_pos, other_paralog_same_pos, other_paralog_other_pos
        if pos == 0:
            return 1, 2, 3
        elif pos == 1:
            return 0, 3, 2
        elif pos == 2:
            return 3, 0, 1
        elif pos == 3:
            return 2, 1, 0
        else:
            sys.exit('Position out of range.')


if __name__ == '__main__':

    pm_model = 'HKY'
    x_js = np.concatenate((np.log([0.3, 0.4, 0.2, 9.5, 1.2, 2.5]), np.log([0.3, 1.0 / 30.0 ])))
    IGC_pm = 'One rate'
    n_orlg = 3
    force = None
    force = {6:0.0}
    rate_variation = True
    test = PSJSModel(x_js, pm_model, n_orlg, IGC_pm, rate_variation, force)
    self = test

    transition = [(0, 2, 2, 3), (0, 2, 2, 1)]
    state_from, state_to = transition
    IGC_init, IGC_p = np.exp(self.x_IGC)
    n = 50


    IGC_0_not_n = IGC_init / IGC_p * (1 - (1 - IGC_p) ** n)
    IGC_0_and_n = IGC_init / IGC_p * (1 - IGC_p) ** n
    print IGC_0_not_n, IGC_0_and_n
    print test.is_transition_compatible(transition), test.cal_IGC_transition_rate(transition, n, (1, 1)), test.cal_IGC_transition_rate(transition, n, (1, 3))
    print test.PMModel.Q_mut
    a = test.get_IGC_process_definition(10)
    print len(a['row_states'])
    b = test.get_PM_process_definition((1, 2))
    c = test.get_PM_process_definition((3, 1))
    print b == c


    
    for i in range(4):
        for j in range(4):
            if i == j:
                continue

            for m in [1,2,3]:
                assert(len(set([test.cal_PM_transition_rate(([k,i,k,i], [k,j,k,j]),(L,m)) for k in range(4) for L in [1,2,3]] + \
                               [test.cal_PM_transition_rate(([i,k,i,k], [j,k,j,k]),(m,L)) for k in range(4) for L in [1,2,3]] )))
                           
            
            print i, j, test.cal_PM_transition_rate(([0,i,0,i], [0,j,0,j]),(1,1))

    
    state_from = (0,1,0,3)
##    for i in range(4):
##        for j in range(4):
##            for k in range(4):
##                for l in range(4):
##                    
##                    state_to = (i,j,k,l)
##                    num_diff = sum([state_from[it] != state_to[it] for it in range(4)])
##                    if num_diff > 3:
##                        continue
##                    elif num_diff == 1:
##                        for m in range(1, 4):
##                            for mm in range(1, 4):
##                                print state_from, state_to, (m, mm), test.cal_IGC_transition_rate([state_from, state_to], n, (m,mm)),\
##                                      test.cal_IGC_transition_rate([state_from, state_to], n, (1,1)), \
##                                      test.cal_IGC_transition_rate([state_from, state_to], n, (1,1)) * 1.2,\
##                                      test.cal_IGC_transition_rate([state_from, state_to], n, (1,1)) * 2.5,\
##                                      (test.cal_IGC_transition_rate([state_from, state_to], n, (1,1)) - IGC_0_not_n)*1.2 + IGC_0_not_n,\
##                                      (test.cal_IGC_transition_rate([state_from, state_to], n, (1,1)) - IGC_0_not_n)*2.5 + IGC_0_not_n
##                        #continue
##                    elif num_diff == 2:
##                        continue
##                        for m in range(1, 4):
##                            for mm in range(1, 4):
##                                if state_to[0] == state_to[2] == state_from[0] and state_to[1] == state_to[3] == state_from[1]:
##                                    print state_from, state_to, test.cal_IGC_transition_rate([state_from, state_to], n, (m,mm))
##                                elif state_to[0] == state_to[2] == state_from[2] and state_to[1] == state_to[3] == state_from[3]:
##                                    print state_from, state_to, test.cal_IGC_transition_rate([state_from, state_to], n, (m,mm))
##                                else:
##                                    continue
##                                    #print state_from, state_to, test.is_transition_compatible([state_from, state_to])
##                    



    
    
