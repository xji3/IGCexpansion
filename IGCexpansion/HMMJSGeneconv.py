# A separate file for joint estimation of HMM + JSIGC model
# Xiang Ji
# xji3@ncsu.edu

from HMMTract import *
from IndCodonGeneconv import IndCodonGeneconv
from copy import deepcopy

class HMMJSGeneconv:
    auto_save_step = 1
    def __init__(self, save_file,
                 # First, ReCodonGeneconv objects
                 newicktree, alignment_file, paralog, summary_path, x, save_path,
                 # Now HMMTract objects
                 IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file,
                 State_List, seq_index_file,
                 model = 'MG94',
                 force = None, nsites = None, clock = None):
        
        self.MG94_IGC = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = force, clock = clock, save_name = save_path + '_'.join(paralog) + '_MG94_IGC_HMMJSGeneconv_save.txt')
        self.MG94_IGC.update_by_x(x[:-1])
        self.MG94_IGC._loglikelihood2()
        
        self.x= x

        self.IGC_sitewise_lnL_file = IGC_sitewise_lnL_file
        self.NOIGC_sitewise_lnL_file = NOIGC_sitewise_lnL_file
        self.MG94_IGC.get_sitewise_loglikelihood_summary(IGC_sitewise_lnL_file, False)
        self.MG94_IGC.get_sitewise_loglikelihood_summary(NOIGC_sitewise_lnL_file, True)
        
        
        self.outgroup_branch = [edge for edge in self.MG94_IGC.edge_list if edge[0] == 'N0' and edge[1] != 'N1'][0]
        Total_blen = sum([self.MG94_IGC.edge_to_blen[edge] for edge in self.MG94_IGC.edge_list if edge != self.outgroup_branch])

        #print 'a'
        self.hmmtract = HMMTract(IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file, State_List, Total_blen, self.MG94_IGC.tau, seq_index_file)
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
        

        # update emission probability
        self.MG94_IGC.get_sitewise_loglikelihood_summary(self.IGC_sitewise_lnL_file, False)
        self.MG94_IGC.get_sitewise_loglikelihood_summary(self.NOIGC_sitewise_lnL_file, True)

        self.hmmtract.IGC_sitewise_lnL = self.hmmtract.read_lnL(self.IGC_sitewise_lnL_file)
        self.hmmtract.NOIGC_sitewise_lnL = self.hmmtract.read_lnL(self.NOIGC_sitewise_lnL_file)

        # update total branch length
        Total_blen = sum([self.MG94_IGC.edge_to_blen[edge] for edge in self.MG94_IGC.edge_list if edge != self.outgroup_branch])
        self.hmmtract.L = Total_blen
        
        # update IGC parameter
        hmm_x = [np.log(self.MG94_IGC.tau) + x[-1], x[-1]]
        self.hmmtract.update_by_x(hmm_x)       

    def _loglikelihood(self, x):
        self.update_by_x(x)
        hmm_x = [x[5] + x[-1], x[-1]]
        ll = self.hmmtract.Forward(False, hmm_x)
        return ll

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
            self.x[:-1] = deepcopy(self.MG94_IGC.x)
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

        

    

        

if __name__ == '__main__':
    pair = ["EDN", "ECP"]
    paralog = pair
    model = 'MG94'
    force = None
    nsites = None
    clock = None

    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
    newicktree = '../test/input_tree.newick'
    Force = None
    output_ctrl = ''
    summary_mat = []
    save_file = '../test/save/HMMJS_' + '_'.join(paralog) + '_MG94_nonclock_save.txt'

    x = np.array([-6.974760482698901809e-01,
                  -5.370097148062192849e-01,
                  -7.240186071588703420e-01,
                  7.278666576839283309e-01,
                  -1.661538873324752974e-01,
                  -6.389150587062580877e-01,
                  -1.627854846878343364e+00,
                  -1.118866867842108759e+00,
                  -3.421887713565998190e+00,
                  -1.856026261560800084e+00,
                  -3.362580076646019211e+00,
                  -2.407034335517133528e+00,
                  -4.244048731094335558e+00,
                  -4.133421357939401020e+00,
                  0.0#-2.6929935812559425#-2.5
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
    NOIGC_sitewise_lnL_file = '../test/Summary/NOIGC_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
    save_path = '../test/save/'

    seq_index_file = '../test/' + '_'.join(paralog) + '_seq_index.txt'

    test = HMMJSGeneconv(save_file, newicktree, alignment_file, paralog, summary_path, x, save_path, IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file,
                         state_list, seq_index_file)

    self = test
    print test._loglikelihood(x)
    test.get_mle(display = True, two_step = False)

    
