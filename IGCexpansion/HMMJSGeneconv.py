# A separate file for joint estimation of HMM + JSIGC model
# Xiang Ji
# xji3@ncsu.edu
import sys
sys.path.append('/usr/local/lib/python2.7/site-packages')
from HMMTract import *
from IndCodonGeneconv import IndCodonGeneconv
from copy import deepcopy
import numdifftools as nd

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
        # x = np.concatenate((self.MG94_IGC.x, np.log(tract_p)))

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
        self.hmmtract.update_by_x([self.x[5], self.x[-1]])
        

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

    def get_mle(self, display = True, two_step = True, One_Dimension = False):
        if two_step:
            self.MG94_IGC.get_mle()
            self.x[:-1] = deepcopy(self.MG94_IGC.x)
            self.update_by_x(self.x)

        if One_Dimension:
            result = self.hmmtract.get_mle(display, True)
            self.x[5] = self.hmmtract.x[0] - self.hmmtract.x[1]
            self.x[-1] = self.hmmtract.x[-1]
            self.update_by_x(self.x)
        else:
            self._loglikelihood(self.x)
            f = partial(self.objective, display)
            guess_x = deepcopy(self.x)
            bnds = [(None, 0.0)] * 3 + [(None, None)] * (len(self.MG94_IGC.x) -3) + [(None, 0.0)]
            result = scipy.optimize.minimize(f, guess_x, jac = False, method = 'L-BFGS-B', bounds = bnds)
            self.save_x()

        if display:
            print(result)

        return result

    def save_x(self):
        np.savetxt(open(self.save_file, 'w+'), self.x.T)
        self.MG94_IGC.save_x()

    def get_summary(self, summary_file):
        out = [self.hmmtract.Forward(False, self.hmmtract.x)]
        out.extend(self.MG94_IGC.pi)
        out.extend([self.MG94_IGC.kappa, self.MG94_IGC.omega, self.MG94_IGC.tau])
        out.extend(np.exp(self.hmmtract.x))
        label = ['ll', 'pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'omega', 'tau', 'eta', 'tract_p']
        k = len(label)  # record the length of non-blen parameters

        # Now add in branch length summary
        label.extend(self.MG94_IGC.edge_list)
        out.extend([self.MG94_IGC.edge_to_blen[label[j]] for j in range(k, len(label))])

        # Now convert labels
        for i in range(k, len(label)):
            label[i] = '(' + ','.join(label[i]) + ')'
            
        footer = ' '.join(label)
        np.savetxt(open(summary_file, 'w+'), np.matrix(out).T, delimiter = ' ', footer = footer) 

    
    def plot_tract_p(self, log_p_list, plot_file):
        self.MG94_IGC.get_mle()
        self.x[:-1] = deepcopy(self.MG94_IGC.x)
        self.update_by_x(self.x)
        ll_list = []
        for log_p in log_p_list:
            ll = -self.hmmtract.objective_1D(False, [log_p])
            ll_list.append(ll)

        with open(plot_file, 'w+') as f:
            f.write('# log_p \t lnL \t \n')
            for it in range(len(log_p_list)):
                f.write('\t'.join([str(log_p_list[it]), str(ll_list[it]), '\n']))

    def get_Hessian(self, One_Dimension = False):
        if One_Dimension:
            f = nd.Derivative(partial(self.hmmtract.objective_1D, False), n = 2)
            result = -f(self.hmmtract.x[1:])
        else:
            f = nd.Hessian(self._loglikelihood)
            result = -f(self.x)
        return result

        
    

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

    summary_file_1D = '../test/Summary/HMM_' + '_'.join(paralog) + '_MG94_nonclock_1D_summary.txt'
    summary_file_all_Dimension = '../test/Summary/HMM_' + '_'.join(paralog) + '_MG94_nonclock_all_summary.txt'

    log_p_list = np.log(3.0/np.array(range(3, 501)))
    plot_file = '../test/plot/HMM_' + '_'.join(paralog) + '_lnL_1D_surface.txt'
    test.plot_tract_p(log_p_list, plot_file)
    test.get_mle(display = True, two_step = True, One_Dimension = True)
    test.get_summary(summary_file_1D)
##    test.get_mle(display = True, two_step = False, One_Dimension = False)
##    test.get_summary(summary_file_all_Dimension)


    


