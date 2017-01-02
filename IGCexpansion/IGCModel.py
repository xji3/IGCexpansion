# A separate class to represent IGC process model
# Xiang Ji
# xji3@ncsu.edu
import numpy as np

class IGCModel:
    supported = ['One rate', 'Most general', 'Symmetric general']          # supported IGC parameterization
    def __init__(self, x_IGC, n_ortholog, parameterization):
        self.pm             = parameterization      # name of parameterization
        self.x_IGC          = x_IGC                 # an array of log() values
        self.n_orlg         = n_ortholog            # total number of ortholog groups

        self.Q_IGC          = None                  # IGC process matrix of size
        self.parameters     = dict()                # a dictionary used to store all IGC parameters
        self.parameter_list = None

        self.init_Q()

    def init_Q(self):
        assert(self.pm in self.supported)  # check if parameterization is implemented
        if self.pm == 'One rate':
            self.init_one_rate_Q()
        if self.pm == 'Most general':
            self.init_most_general_Q()
        if self.pm == 'Symmetric general':
            self.init_sym_general_Q()

    def init_one_rate_Q(self):
        assert(len(self.x_IGC) == 1)  # check x_IGC length first
        self.Q_IGC = np.ones((self.n_orlg, self.n_orlg), dtype = np.floating) * np.exp(self.x_IGC[0])
        np.fill_diagonal(self.Q_IGC, 0.0)
        self.parameters['Tau'] = np.exp(self.x_IGC[0])

    def init_most_general_Q(self):
        # TODO: finish self.parameters in this case
        assert(len(self.x_IGC) == self.n_orlg * (self.n_orlg - 1))
        # x_IGC in this case should be in linear order of
        # t_{1,2}, t_{1, 3}, ..., t_{1, n}, t_{2, 1}, ... t_{2, n}, ..., t{n, n-1}
        IGC_array = []
        IGC_rates = np.exp(self.x_IGC)
        for n_it in range(self.n_orlg):
            row = IGC_rates[n_it * (self.n_orlg - 1) : (n_it + 1) * (self.n_orlg - 1)]            
            IGC_array.append(np.insert(row, n_it, 0.0))
        self.Q_IGC = np.array(IGC_array)

    def init_sym_general_Q(self):
        assert(len(self.x_IGC) == self.n_orlg * (self.n_orlg - 1) / 2)
        # x_IGC in this case should be in linear order of
        # t_{1,2}, t_{1, 3}, ..., t_{1, n}, t_{2, 3}, ... t_{2, n}, ..., t{n - 1, n}
        IGC_array = []

    def __str__(self): # overide for print function
        return 'IGC Model parameterization: ' + self.pm + '\n' + 'IGC parameters: ' + str(np.exp(self.x_IGC))
            

if __name__ == '__main__':
    test = IGCModel(np.log([4.9]), 7, 'One rate')
    print test.Q_IGC
    test = IGCModel(np.log(range(1, 7)), 3, 'Most general')
    print test.Q_IGC
    self = test
