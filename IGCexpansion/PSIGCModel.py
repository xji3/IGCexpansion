# A separate class to represent IGC process model
# Xiang Ji
# xji3@ncsu.edu
import numpy as np

class PSIGCModel:
    supported = ['One rate']          # supported IGC parameterization
    def __init__(self, x_IGC, n_ortholog, parameterization, accessible_orlg_pair = None, force = None):
        self.pm             = parameterization      # name of parameterization
        self.x_IGC          = x_IGC                 # an array of log() values
        self.n_orlg         = n_ortholog            # total number of ortholog groups
        self.force          = force                 # used for parameter value constraint

        self.rate_IGC       = None                  # IGC initiation rate matrix
        self.p_IGC          = None                  # IGC tract length parameter p matrix
        
        self.parameters     = dict()                # a dictionary used to store all IGC parameters
        self.parameter_list = list()

        self.accessible_orlg_pair = accessible_orlg_pair   # list of (row, state) indices of accessible orlg_pair in Q_IGC, used for assign/update Q_IGC

        self.init_Q()

    def init_Q(self):
        assert(self.pm in self.supported)  # check if parameterization is implemented
        if self.pm == 'One rate':
            self.init_one_rate_Q()
        else:
            sys.exit('Parameterization not implemented in PSIGCModel class!')


    def update_by_x_IGC(self, new_x_IGC):
        assert(len(self.x_IGC) == len(new_x_IGC))
        self.x_IGC = new_x_IGC
        self.init_Q()
    
    def init_one_rate_Q(self):
        assert(len(self.x_IGC) == 2 and self.x_IGC[1] <= 0.0)
        # check x_IGC length first and the tract length geometric distribution parameter p should be <= 1.0
        if not self.force == None:
            assert(all([ -1 < key < 2 for key in self.force]))

        IGC_rate = np.exp(self.x_IGC)
        if not self.force == None: # add in force constraints
            IGC_rate = self.force[0]
        init_rate, tract_p = IGC_rate
        
        self.rate_IGC = np.ones((self.n_orlg, self.n_orlg), dtype = np.floating) * init_rate
        self.p_IGC = np.ones((self.n_orlg, self.n_orlg), dtype = np.floating) * tract_p
        
        np.fill_diagonal(self.rate_IGC, 0.0)
        np.fill_diagonal(self.p_IGC, 0.0)
        
        self.parameters['Init_rate'] = init_rate
        self.parameters['tract_p'] = tract_p
        if len(self.parameter_list) != 2:
            self.parameter_list = ['Init_rate', 'tract_p']

               


    def __str__(self): # overide for print function
        return 'IGC Model parameterization: ' + self.pm + '\n' + 'IGC parameters: ' + str(np.exp(self.x_IGC)) \
               + '\n Accessible orlg pairs:' + str(self.accessible_orlg_pair)
            

if __name__ == '__main__':
    test = PSIGCModel(np.log([4.9, 0.5]), 7, 'One rate', {0:0.0})
    print test.rate_IGC
    print test.p_IGC

