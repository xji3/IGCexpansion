# A separate class to represent IGC process model
# Xiang Ji
# xji3@ncsu.edu
import numpy as np

class IGCModel:
    supported = ['One rate', 'Most general', 'Symmetric general']          # supported IGC parameterization
    def __init__(self, x_IGC, n_ortholog, parameterization, accessible_orlg_pair = None, force = None):
        self.pm             = parameterization      # name of parameterization
        self.x_IGC          = x_IGC                 # an array of log() values
        self.n_orlg         = n_ortholog            # total number of ortholog groups
        self.force          = force                 # used for parameter value constraint

        self.Q_IGC          = None                  # IGC process matrix of size
        self.parameters     = dict()                # a dictionary used to store all IGC parameters
        self.parameter_list = list()

        self.accessible_orlg_pair = accessible_orlg_pair   # list of (row, state) indices of accessible orlg_pair in Q_IGC, used for assign/update Q_IGC

        self.init_Q()

    def init_Q(self):
        assert(self.pm in self.supported)  # check if parameterization is implemented
        if self.pm == 'One rate':
            self.init_one_rate_Q()
        if self.pm == 'Most general':
            self.init_most_general_Q()
        if self.pm == 'Symmetric general':
            self.init_sym_general_Q()

    def update_by_x_IGC(self, new_x_IGC):
        assert(len(self.x_IGC) == len(new_x_IGC))
        self.x_IGC = new_x_IGC
        self.init_Q()
    
    def init_one_rate_Q(self):
        assert(len(self.x_IGC) == 1)  # check x_IGC length first
        if not self.force == None:
            assert(all([ -1 < key < 1 for key in self.force]))

        IGC_rate = np.exp(self.x_IGC[0])
        if not self.force == None: # add in force constraints
            IGC_rate = self.force[0]
        self.Q_IGC = np.ones((self.n_orlg, self.n_orlg), dtype = np.floating) * IGC_rate
        np.fill_diagonal(self.Q_IGC, 0.0)
        self.parameters['Tau'] = IGC_rate
        self.parameter_list = ['Tau']

    def init_most_general_Q(self):
        assert(self.accessible_orlg_pair)
        assert(len(self.x_IGC) == 2 * len(self.accessible_orlg_pair) <= self.n_orlg * (self.n_orlg - 1) )
        if not self.force == None:
            assert(all([-1 < key < len(self.x_IGC) for key in self.force]))
        # x_IGC in this case should be in linear order of
        # t_{1,2}, t_{1, 3}, ..., t_{1, n}, t_{2, 1}, ... t_{2, n}, ..., t{n, n-1}
        IGC_array = []
        IGC_rates = np.exp(self.x_IGC)

        # Now update force constraints
        if not self.force == None:
            for i in self.force.keys():
                IGC_rates[i] = self.force[i]

        self.Q_IGC = np.zeros((self.n_orlg, self.n_orlg), dtype = float)
        for idx in range(len(self.accessible_orlg_pair)):
            i, j = self.accessible_orlg_pair[idx]
            self.Q_IGC[i, j] = IGC_rates[idx]
            self.Q_IGC[j, i] = IGC_rates[idx + len(self.accessible_orlg_pair)]


        # Now push the rates into self.parameters
        for idx in range(len(self.accessible_orlg_pair)):
            i, j = self.accessible_orlg_pair[idx]
            if not 'Tau_' + str(i) + '->' + str(j) in self.parameter_list:
                self.parameter_list.append('Tau_' + str(i) + '->' + str(j))
            self.parameters['Tau_' + str(i) + '->' + str(j)] = self.Q_IGC[i, j]
            if not 'Tau_' + str(j) + '->' + str(i) in self.parameter_list:
                self.parameter_list.append('Tau_' + str(j) + '->' + str(i))
            self.parameters['Tau_' + str(j) + '->' + str(i)] = self.Q_IGC[j, i]

            
    def init_sym_general_Q(self):
        assert(self.accessible_orlg_pair)
        assert(len(self.x_IGC) == len(self.accessible_orlg_pair) <= self.n_orlg * (self.n_orlg - 1) / 2)
        if not self.force == None:
            assert(all([-1 < key < len(self.x_IGC) for key in self.force]))
        # x_IGC in this case should be in linear order of
        # t_{1,2}, t_{1, 3}, ..., t_{1, n}, t_{2, 3}, ... t_{2, n}, ..., t{n - 1, n}
        IGC_array = []
        IGC_rates = np.exp(self.x_IGC)

        # Now update force constraints
        if not self.force == None:
            for i in self.force.keys():
                IGC_rates[i] = self.force[i]

        self.Q_IGC = np.zeros((self.n_orlg, self.n_orlg), dtype = float)
        for idx in range(len(self.accessible_orlg_pair)):
            i, j = self.accessible_orlg_pair[idx]
            self.Q_IGC[i, j] = self.Q_IGC[j, i] = IGC_rates[idx]

        # Now push the rates into self.parameters
        for idx in range(len(self.accessible_orlg_pair)):
            i, j = self.accessible_orlg_pair[idx]
            if not 'Tau_' + str(i) + '<->' + str(j) in self.parameter_list:
                self.parameter_list.append('Tau_' + str(i) + '<->' + str(j))
            self.parameters['Tau_' + str(i) + '<->' + str(j)] = self.Q_IGC[i, j]
                
                


    def __str__(self): # overide for print function
        return 'IGC Model parameterization: ' + self.pm + '\n' + 'IGC parameters: ' + str(np.exp(self.x_IGC)) \
               + '\nAccessible orlg pairs:' + str(self.accessible_orlg_pair)
            

if __name__ == '__main__':
    test = IGCModel(np.log([4.9]), 7, 'One rate', {0:0.0})
    print test.Q_IGC


    from Tree import Tree
    from Common import *
    from Func import *

    tree_newick = '../test/PrimateTest.newick'
    DupLosList = '../test/PrimateTestDupLost.txt'


    node_to_pos = {'D1':0, 'D2':0, 'D3':1, 'D4':2, 'L1':2}
    terminal_node_list = ['Chinese_Tree_Shrew', 'Macaque', 'Olive_Baboon', 'Orangutan', 'Gorilla', 'Human']
    tree = Tree(tree_newick, DupLosList, terminal_node_list, node_to_pos)

    conf_list = count_process(tree.node_to_conf)
    accessible_orlg_pair = get_accessible_orlg_pair(conf_list)

    test = IGCModel(np.log(range(1, 1 + len(accessible_orlg_pair))), tree.n_orlg, 'Symmetric general', accessible_orlg_pair)
    print test.Q_IGC
    self = test

    x_IGC = np.log(range(2, 2 + len(accessible_orlg_pair) * 2))
    test = IGCModel(x_IGC, tree.n_orlg, 'Most general', accessible_orlg_pair)
    print test.Q_IGC
    self = test
