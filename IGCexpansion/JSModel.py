# A separate class to represent Joint State IGC model (JS IGC models)
# JS IGC model = IGC model + Point mutation model
# Xiang Ji
# xji3@ncsu.edu

from IGCModel import IGCModel
from PMModel import PMModel
import numpy as np
import itertools
from copy import deepcopy

class JSModel:
    def __init__(self, n_js, x_js, pm_model, n_orlg, IGC_pm):
        self.n_js   = n_js            # number of contemporaneous paralog states considered on each branch
        self.x_js   = x_js            # one concatenated vector to store all rate matrix parameters
        self.x_pm   = None            # x_pm vector for PMModel
        self.x_IGC  = None            # x_IGC vector for IGCModel

        self.pm_model = pm_model      # name of point mutation model
        self.IGC_pm   = IGC_pm        # IGC parameterization
        self.n_orlg   = n_orlg        # total number of ortholog groups

        self.PMModel  = None          # PMModel class instance for point mutation model
        self.IGCModel = None          # IGCModel class instance for IGC model

        self.state_space_shape = None # initialized in init_models() function

        self.init_models()


    def unpack_x_js(self):
        # first, check if the models are supported
        assert(self.pm_model in PMModel.supported)
        assert(self.IGC_pm in IGCModel.supported)
        if self.pm_model == 'HKY':
            num_x_pm = 4
        else:
            print 'The point mutation model is not supported.'
            num_x_pm = 0

        if self.IGC_pm == 'One rate':
            num_x_IGC = 1
        elif self.IGC_pm == 'Most general':
            num_x_IGC = self.n_orlg * (self.n_orlg - 1)
        else:
            print 'The IGC parameterization has not been implemented.'
            num_x_IGC = 0

        self.x_pm = self.x_js[:num_x_pm]
        self.x_IGC = self.x_js[num_x_pm:]
        assert(num_x_pm + num_x_IGC == len(self.x_js))

    def init_models(self):
        self.unpack_x_js()
        if self.pm_model == 'HKY':
            self.state_space_shape = [4 for i in range(self.n_js)]
        else:
            print 'The point mutation model has not been implemented.'
        self.PMModel = PMModel(self.pm_model, self.x_pm)
        self.IGCModel = IGCModel(self.x_IGC, self.n_orlg, self.IGC_pm)
        print(self.PMModel)
        print(self.IGCModel)
        
    def is_configuration(self, configuration):
        return len(configuration) == self.n_js and\
               all([len(single_conf) == 2 for single_conf in configuration]) and \
               all([single_conf[1] == 0 or single_conf[1] == 1 for single_conf in configuration]) and\
               all([ -1 < configuration[i][0] < self.n_orlg for i in range(self.n_js)])

    def divide_configuration(self, configuration):
        assert(self.is_configuration(configuration))
        ortho_group_to_pos = dict(extent = {}, distinct = [])
        # extent positions that represent same paralog (by the same ortholog group number) have to be in the same state
        # distinct positions don't change states, thus only track positions
        for pos in range(self.n_js):
            if configuration[pos][1] == 1: # extent
                ortho_group = configuration[pos][0]
                if ortho_group in ortho_group_to_pos['extent']:
                    ortho_group_to_pos['extent'][ortho_group].append(pos)
                else:
                    ortho_group_to_pos['extent'][ortho_group] = [pos]
            elif configuration[pos][1] == 0: # distinct
                ortho_group_to_pos['distinct'].append(pos)

        return ortho_group_to_pos
            
            
    def is_state_compatible(self, state, configuration):
        assert(self.is_configuration(configuration)) # it has to be a right configuration first

        # The length of state should be self.n_js
        indicator = (len(state) == self.n_js)

        ortho_group_to_pos = self.divide_configuration(configuration)
        for ortho_group in ortho_group_to_pos['extent'].keys():
            pos_list = ortho_group_to_pos['extent'][ortho_group]
            indicator = indicator and len(set([state[pos] for pos in pos_list])) == 1 \
                        and len(set([self.state_space_shape[pos] for pos in pos_list])) == 1

        for pos in range(self.n_js):
            indicator = indicator and -1 < state[pos] < self.state_space_shape[pos]
        return indicator

    def is_transition_compatible(self, transition, configuration):
        assert(len(transition) == 2) # transition should contain two states
        state_from, state_to = transition
        if state_from == state_to:
            return False
        indicator = self.is_state_compatible(state_from, configuration)\
                    and self.is_state_compatible(state_to, configuration)

        if not indicator:
            return False
        # Two states should only differ at one ortholog group
        pos_list = [i for i in range(self.n_js) if state_from[i] != state_to[i]]
        indicator = indicator and len(set([state_from[i] for i in pos_list])) == 1\
                    and len(set([configuration[i][0] for i in pos_list if configuration[i][1] == 1])) == 1\
                    and state_from[pos_list[0]] != state_to[pos_list[0]]

        if not indicator:
            return False
        # Distinct lineage should not change state
        ortho_group_to_pos = self.divide_configuration(configuration)
        for pos in ortho_group_to_pos['distinct']:
            indicator = indicator and state_from[pos] == state_to[pos]
        return indicator

    def cal_js_transition_rate(self, transition, configuration):
        # this function should only be called for allowed transitions
        # this step may be redundant, but it doesnot cost much
        assert(self.is_transition_compatible(transition, configuration))

        # now need to find the position
        state_from, state_to = transition
        pos_list = [i for i in range(self.n_js) if state_from[i] != state_to[i]]

        # get transition rate from Point Mutation Model
        q_ij = self.PMModel.Q_mut[state_from[pos_list[0]], state_to[pos_list[0]]]

        # Now get IGC rates from IGC Model
        ortho_group_to_pos = self.divide_configuration(configuration)
        IGC_ortho_groups = [ortho_group for ortho_group in ortho_group_to_pos['extent'] \
                            if state_from[ortho_group_to_pos['extent'][ortho_group][0]] == \
                            state_to[pos_list[0]]]
        t_IGC = 0.0
        IGC_to = configuration[pos_list[0]][0]
        if len(IGC_ortho_groups) > 0:
            for ortho_group in IGC_ortho_groups:
                t_IGC += self.IGCModel.Q_IGC[ortho_group, IGC_to]

        return q_ij + t_IGC

    def get_js_transition_rates_BF(self, configuration):  # Brute force way for test purpose
        for state_from in itertools.product(range(self.state_space_shape[0]), repeat = self.n_js):
            for state_to in itertools.product(range(self.state_space_shape[0]), repeat = self.n_js):
                if self.is_transition_compatible([state_from, state_to], configuration):
                    yield state_from, state_to, self.cal_js_transition_rate([state_from, state_to], configuration)

    def get_js_transition_rates(self, configuration):
        ortho_group_to_pos = self.divide_configuration(configuration)
        for state_from in itertools.product(range(self.state_space_shape[0]), repeat = self.n_js):
            for orlg in ortho_group_to_pos['extent']:
                state_to = list(deepcopy(state_from))
                for nt in range(self.state_space_shape[ortho_group_to_pos['extent'][orlg][0]]):
                    for pos in ortho_group_to_pos['extent'][orlg]:
                        state_to[pos] = nt                            
                        if self.is_transition_compatible([state_from, state_to], configuration):
                            yield state_from, state_to, self.cal_js_transition_rate([state_from, state_to], configuration)


    def get_process_definition(self, configuration):
        row_states = []
        column_states = []
        transition_rates = []
        for row_state, col_state, transition_rate in self.get_js_transition_rates(configuration):
            row_states.append(row_state)
            column_states.append(col_state)
            transition_rates.append(transition_rate)
        process_definition = dict(
            row_states = row_states,
            column_states = column_states,
            transition_rates = transition_rates)
        return process_definition

        
        


if __name__ == '__main__':
    pm_model = 'HKY'
    x_js = np.log([0.3, 0.5, 0.2, 9.5, 4.9])
    n_orlg = 4
    IGC_pm = 'One rate'
    n_js = 5
    test = JSModel(n_js, x_js, pm_model, n_orlg, IGC_pm)
    self = test

    test_configuration = [(i/2, 1) for i in range(n_js)]
    print 'test configuration = ', test_configuration
    print test.is_configuration(test_configuration)
    print test.divide_configuration(test_configuration)
    print test.is_state_compatible([2, 2, 3, 3, 1], test_configuration),\
          test.is_state_compatible([2, 2, 2, 3, 1], test_configuration),\
          test.is_state_compatible([2, 2, 4, 4, 1], test_configuration)

    transition = ([2, 2, 3, 3, 1], [2, 2, 4, 4, 1])
    print test.is_transition_compatible(transition, test_configuration)
    transition = ([2, 2, 3, 3, 1], [2, 2, 0, 0, 1])
    print test.cal_js_transition_rate(transition, test_configuration)
    transition = ([2, 2, 3, 3, 1], [2, 2, 2, 2, 1])
    print test.cal_js_transition_rate(transition, test_configuration)

    test_configuration_2p = [(0, 1), (0, 1), (1, 1), (1, 1), (0, 1)]
    process_definition = test.get_process_definition(test_configuration)
    print len(process_definition['row_states'])
    #transition = [[0, 0, 0, 0, 0], [1, 1, 1, 1, 0]]
    #test.cal_js_transition_rate(transition, test_configuration_2p)

    
