# A separate file for joint estimation of HMM + JSIGC model
# Xiang Ji
# xji3@ncsu.edu

from HMMTract import *
from CodonGeneconv import ReCodonGeneconv
from copy import deepcopy

class HMMJSGeneconv:
    auto_save_step = 1
    def __init__(self, save_file,
                 # First, ReCodonGeneconv objects
                 newicktree, alignment_file, paralog, summary_path, x, save_path,
                 # Now HMMTract objects
                 IGC_sitewise_lnL_file, Force_sitewise_lnL_file,
                 State_List, seq_index_file,
                 model = 'MG94',
                 force = None, nsites = None, clock = None):
        
        self.MG94_IGC = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = force, clock = clock, save_path = save_path)
        self.MG94_Force = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = {5:0.0}, clock = clock, save_path = save_path)
        self.MG94_IGC.update_by_x(x[:-1])
        self.MG94_Force.update_by_x(x[:-1])
        self.MG94_IGC._loglikelihood2()
        self.MG94_Force._loglikelihood2()
        
        self.x= x

        self.IGC_sitewise_lnL_file = IGC_sitewise_lnL_file
        self.Force_sitewise_lnL_file = Force_sitewise_lnL_file
        self.MG94_IGC.get_sitewise_loglikelihood_summary(IGC_sitewise_lnL_file)
        self.MG94_Force.get_sitewise_loglikelihood_summary(Force_sitewise_lnL_file)
        
        self.outgroup_branch = [edge for edge in self.MG94_IGC.edge_list if edge[0] == 'N0' and edge[1] != 'N1'][0]
        Total_blen = sum([self.MG94_IGC.edge_to_blen[edge] for edge in self.MG94_IGC.edge_list if edge != self.outgroup_branch])
        
        self.hmmtract = HMMTract(IGC_sitewise_lnL_file, Force_sitewise_lnL_file, State_List, Total_blen, self.MG94_IGC.tau, seq_index_file)
        self.hmmtract.update_by_x(np.log([self.MG94_IGC.tau * 0.05, 0.05]))

        #self.update_by_x(self.x)
        self.save_file = save_file
        if os.path.isfile(self.save_file):
            self.initialize_by_save()
            print 'Loaded parameter values from save file: ', self.save_file

        self.auto_save = 0

    def initialize_by_save(self):
        self.x = np.loadtxt(open(self.save_file, 'r'))
        self.update_by_x(self.x)

    def update_by_x(self, x):
        self.x = deepcopy(x)
        # update two MG94 + IGC models first
        self.MG94_IGC.update_by_x(x[:-1])
        self.MG94_Force.update_by_x(x[:-1])

        # update emission probability
        self.MG94_IGC.get_sitewise_loglikelihood_summary(self.IGC_sitewise_lnL_file)
        self.MG94_Force.get_sitewise_loglikelihood_summary(self.Force_sitewise_lnL_file)
        
        self.hmmtract.IGC_sitewise_lnL = self.hmmtract.read_lnL(self.IGC_sitewise_lnL_file)
        self.hmmtract.Force_sitewise_lnL = self.hmmtract.read_lnL(self.Force_sitewise_lnL_file)

        # update total branch length
        Total_blen = sum([self.MG94_IGC.edge_to_blen[edge] for edge in self.MG94_IGC.edge_list if edge != self.outgroup_branch])
        self.hmmtract.L = Total_blen
        
        # update IGC parameter
        hmm_x = [np.log(self.MG94_IGC.tau) + x[-1], x[-1]]
        self.hmmtract.update_by_x(hmm_x)       

    def _loglikelihood(self, x):
        self.update_by_x(x)
        hmm_x = [x[5] + x[-1], x[-1]]
        ll = -self.hmmtract.Forward(False, hmm_x)
        return -ll

    def objective(self, display, x):
        assert(len(x) == len(self.x))
        if display:
            print 'New x array: ', x
        ll = self._loglikelihood(x)
        if display:
            print 'Current exp(x) array: ', np.exp(self.x)
            print 'lnL = ', -ll

        if self.auto_save == HMMJSGeneconv.auto_save_step * len(self.x):
            self.save_x()
            self.auto_save == 0
        else:
            self.auto_save += 1
            
        return ll

    def get_mle(self, display = True, two_step = True):
        if two_step:
            self.MG94_IGC.get_mle()
            self.MG94_Force.update_by_x(self.MG94_IGC.x)
            self.x[:-1] = deepcoy(self.MG94_IGC.x)
        self._loglikelihood(self.x)
        f = partial(self.objective, display)
        guess_x = deepcopy(self.x)
        bnds = [(None, 0.0)] * 3 + [(None, None)] * (len(self.MG94_IGC.x) -3) + [(None, 0.0)]
        result = scipy.optimize.minimize(f, guess_x, jac = False, method = 'L-BFGS-B', bounds = bnds)

        if display:
            print(result)

    def save_x(self):
        np.savetxt(open(self.save_file, 'w+'), self.x.T)
        self.MG94_IGC.save_x()
        self.MG94_Force.save_x()
        

    

        

if __name__ == '__main__':
    pair = ["EDN", "ECP"]
    paralog = pair

    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
    newicktree = '../test/input_tree.newick'
    Force = None
    output_ctrl = ''
    summary_mat = []
    save_file = '../test/save/HMMJS_' + '_'.join(paralog) + '_MG94_nonclock_save.txt'

    x = np.array([-6.972787845489714087e-01,
                  -5.371080077441253708e-01,
                  -7.240047380894862883e-01,
                  7.238578797129271436e-01,
                  -1.898373482480217311e-01,
                  -1.4183244411092852, #-4.600735470664583104e-01,
                  -1.613942190237640295e+00,
                  -1.123019321880592836e+00,
                  -3.443903685022769334e+00,
                  -1.854709904478881510e+00,
                  -3.365339758069785692e+00,
                  -2.410226233318178313e+00,
                  -4.249147091385517605e+00,
                  -4.137233772266141862e+00,
                  -2.6929935812559425#-2.5
                  ])
    x_2 = np.array([-0.69727878,
                  -0.53710801,
                  -0.72400474,
                  0.72385788,
                  -0.18983735,
                  -2.28199719,
                  -2.01452935,
                  -1.0337633,
                  -3.29369029,
                  -1.75318807,
                  -3.25869777,
                  -2.27341043,
                  -4.20160402,
                  -4.110472,
                  -2.83333001])
    print 
    print '**' + '_'.join(paralog)+ '**', output_ctrl
    
    alignment_file = '../test/EDN_ECP_Cleaned.fasta'
    summary_path = '../test/Summary/'
    IGC_sitewise_lnL_file = '../test/Summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
    Force_sitewise_lnL_file = '../test/Summary/Force_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
    save_path = '../test/save/'

    seq_index_file = '../test/' + '_'.join(paralog) + '_seq_index.txt'

    test = HMMJSGeneconv(save_file, newicktree, alignment_file, paralog, summary_path, x_2, save_path, IGC_sitewise_lnL_file, Force_sitewise_lnL_file,
                         state_list, seq_index_file)

    self = test
    #print test._loglikelihood(x)
    test.get_mle(display = True, two_step = False)

    #HKY_IGC = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = None, clock = None, save_path = save_path)
    #[a, b] = HKY_IGC.get_HKYGeneconv()

