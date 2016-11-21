# A separate class to represent Joint State IGC model (JS IGC models)
# JS IGC model = IGC model + Point mutation model
# Xiang Ji
# xji3@ncsu.edu

from IGCModel import IGCModel
from PMModel import PMModel
import numpy as np

class JSModel:
    def __init__(self, n_js, x_js, pm_model, n_orlg, IGC_pm):
        self.n_js   = n_js          # number of contemporaneous paralog states considered on each branch
        self.x_js   = x_js          # one concatenated vector to store all rate matrix parameters
        self.x_pm   = None          # x_pm vector for PMModel
        self.x_IGC  = None          # x_IGC vector for IGCModel

        self.pm_model = pm_model    # name of point mutation model
        self.IGC_pm   = IGC_pm      # IGC parameterization
        self.n_orlg   = n_orlg      # total number of ortholog groups

        self.PMModel  = None        # PMModel class instance for point mutation model
        self.IGCModel = None        # IGCModel class instance for IGC model

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
        self.PMModel = PMModel(self.pm_model, self.x_pm)
        self.IGCModel = IGCModel(self.x_IGC, self.n_orlg, self.IGC_pm)
        print(self.PMModel)
        print(self.IGCModel)
        
    def is_configuration(self, configuration):
        return len(configuration) == self.n_js and\
               all([len(single_conf) == 2 for single_conf in configuration]) and \
               all([single_conf[1] == 0 or single_conf[1] == 1 for single_conf in configuration])

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
        indicator = True
        ortho_group_to_pos = self.divide_configuration(configuration)
        for ortho_group in ortho_group_to_pos['extent'].keys():
            pos_list = ortho_group_to_pos['extent'][ortho_group]
            indicator = indicator and len(set([state[pos] for pos in pos_list])) == 1
        return indicator
        


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
    print test.is_state_compatible([2, 2, 3, 3, 1, 1], test_configuration),\
          test.is_state_compatible([2, 2, 2, 3, 1, 1], test_configuration)
