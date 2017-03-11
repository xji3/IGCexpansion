# A separate file for IGC model with tract length and initiation rate
# This class is mainly used for simulations
# Xiang Ji
# xji3@ncsu.edu

import numpy as np

class IGCTractModel:
    # initiation rate matrix parameterization
    supported_init_pm = ['One rate', 'Most general', 'Symmetric general']
    # tract length matrix parameterization 
    supported_tract_pm = ['One rate']

    def __init__(self, x_IGC, orlg_list, parameterization):
        self.pm             = parameterization         # list of two names of parameterization
        self.x_IGC          = x_IGC                    # an array of REAL values
        self.orlg_list      = orlg_list                # list of contemporaneous paralogs (their orthologous group numbers)

        self.x_init         = None
        self.x_tract        = None
        self.n_orlg         = len(self.orlg_list)
        self.Q_tract        = None                     # tract length matrix
        self.Q_init         = None                     # initiation rate matrix
        self.parameters     = dict()
        self.parameter_list = list()
        self.num_to_orlg    = {i:sorted(self.orlg_list, reverse=False)[i] for i in range(self.n_orlg)}
        self.orlg_to_num    = {self.num_to_orlg[i]:i for i in range(self.n_orlg)}
        self.init_Q()

    def unpack_x_IGC(self):
        assert(len(self.pm) == 2)
        assert(self.pm[0] in self.supported_init_pm and self.pm[1] in self.supported_tract_pm)
        if self.pm[0] == 'One rate':
            num_init = 1

        if self.pm[1] == 'One rate':
            num_tract = 1

        assert(len(self.x_IGC) == num_init + num_tract)
        self.x_init = self.x_IGC[:num_init]
        self.x_tract = self.x_IGC[num_init:]

    def update_orlg_list(self, new_orlg_list):
        self.orlg_list = new_orlg_list
        self.n_orlg = len(new_orlg_list)
        self.num_to_orlg    = {i:sorted(self.orlg_list, reverse=False)[i] for i in range(self.n_orlg)}
        self.orlg_to_num    = {self.num_to_orlg[i]:i for i in range(self.n_orlg)}
        self.init_Q()

        
    def init_Q(self):
        self.unpack_x_IGC()

        if self.pm[0] == 'One rate':
            self.init_one_rate_init_Q()

        if self.pm[1] == 'One rate':
            self.init_one_rate_tract_Q()

    def init_one_rate_init_Q(self):
        assert(len(self.x_init) == 1)
        self.Q_init = np.ones((self.n_orlg, self.n_orlg), dtype = np.floating) * self.x_init[0]
        np.fill_diagonal(self.Q_init, 0.0)
        self.parameters['Init_rate'] = self.x_init[0]
        if not 'Init_rate' in self.parameter_list:
            self.parameter_list.append('Init_rate')

    def init_one_rate_tract_Q(self):
        assert(len(self.x_tract) == 1)
        self.Q_tract = np.ones((self.n_orlg, self.n_orlg), dtype = np.floating) * self.x_tract[0]
        np.fill_diagonal(self.Q_tract, 0.0)
        self.parameters['Tract_p'] = self.x_tract[0]
        if not 'Tract_p' in self.parameter_list:
            self.parameter_list.append('Tract_p')

    def __str__(self): # overide for print function
        return 'IGC initiation parameterization: ' + self.pm[0] + '\n' + \
               'IGC tract parameterization: ' + self.pm[1] + '\n' + \
               'IGC parameters: ' + ' '.join([item + ' '+ str(self.parameters[item]) for item in self.parameter_list])\
               + '\n'

if __name__ == '__main__':
    init_pm = 'One rate'
    tract_pm = 'One rate'
    pm = [init_pm, tract_pm]
    x_IGC = [2.0, 0.3]
    orlg_list = [1, 2, 5]

    test = IGCTractModel(x_IGC, orlg_list, pm)
    self = test
    
