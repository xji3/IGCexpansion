# A separate class to represent Point mutation model
# Xiang Ji
# xji3@ncsu.edu
import numpy as np


class PMModel:
    supported = ['HKY']                         # supported models
    def __init__(self, model_name, x_pm):
        self.name         = model_name          # name of supported models
        self.x_pm         = x_pm                # an array of log() values

        self.Q_mut        = None                # Point mutation Q matrix
        self.parameters   = dict()              # a dictionary to store all model parameters
        self.init_Q()


    def init_Q(self):
        assert(self.name in self.supported) # check if PM model is implemented
        if self.name == 'HKY':
            self.init_HKY_Q()

    def init_HKY_Q(self):
        # This function initialize Q matrix for HKY85 point mutation model
        # first, check x_pm size
        # x_pm = np.log([%AG, %A, %C, kappa])
        # HKY model has 3 parameters for nt frequency, and 1 parameter kappa for transition/transversion ratio
        assert(len(self.x_pm) == 4)
        self.unpack_frequency()
        kappa = np.exp(self.x_pm[3])
        self.parameters['kappa'] = kappa

        Qbasic = np.array([
            [0, 1.0, kappa, 1.0],
            [1.0, 0, 1.0, kappa],
            [kappa, 1.0, 0, 1.0],
            [1.0, kappa, 1.0, 0],
            ]) * np.array(self.parameters['pi'])
        # assume PM at stationary
        stationary_distn = [ self.parameters['pi']['ACGT'.index(nt)] for nt in 'ACGT' ]
        stationary_distn = np.array(stationary_distn) / sum(stationary_distn)
        expected_rate = np.dot(stationary_distn, Qbasic.sum(axis = 1))
        self.Q_mut = Qbasic / expected_rate

               

    def unpack_frequency(self):
        x_process = np.exp(self.x_pm[:3])
        pi_a = x_process[0] * x_process[1]
        pi_c = (1 - x_process[0]) * x_process[2]
        pi_g = x_process[0] * (1 - x_process[1])
        pi_t = (1 - x_process[0]) * (1 - x_process[2])
        self.parameters['pi'] = [pi_a, pi_c, pi_g, pi_t]        

    def update_by_x_pm(self, new_x_pm):
        assert(len(self.x_pm) == len(new_x_pm))
        self.x_pm = new_x_pm
        self.init_Q()

    def __str__(self): # overide for print function
        return 'Point mutation model: ' + self.name + '\n' + \
               'Point mutation parameters: ' + ' '.join([item + ' '+ str(self.parameters[item]) for item in self.parameters])

if __name__ == '__main__':
    test = PMModel('HKY', np.log([0.3, 0.5, 0.2, 9.5]))
    self = test
    print test.Q_mut
    test.update_by_x_pm(np.log([0.1, 0.9, 0.3, 11.0]))
    print test.Q_mut
