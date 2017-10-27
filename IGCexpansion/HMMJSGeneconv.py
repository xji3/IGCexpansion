# A separate file for joint estimation of HMM + JSIGC model
# Xiang Ji
# xji3@ncsu.edu
import sys
sys.path.append('/usr/local/lib/python2.7/site-packages')
from HMMTract import *
from IndCodonGeneconv import IndCodonGeneconv
from JSGeneconv import JSGeneconv
from copy import deepcopy
import numdifftools as nd

class HMMJSGeneconv:
    auto_save_step = 1
    def __init__(self, save_file,
                 # First, IndCodonGeneconv objects
                 newicktree, alignment_file, paralog, summary_path, x, save_path,
                 # Now HMMTract objects
                 IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file,
                 State_List, seq_index_file,
                 model = 'MG94',
                 force = None, nsites = None, clock = None,
                 rate_variation = False):

        if rate_variation:
            self.MG94_IGC = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = force, clock = clock,\
                                              save_name = save_path + '_'.join(paralog) + '_' + model + '_rv_IGC_HMMJSGeneconv_save.txt',\
                                              rate_variation = True)
        else:
            self.MG94_IGC = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = force, clock = clock,\
                                              save_name = save_path + '_'.join(paralog) + '_' + model + '_IGC_HMMJSGeneconv_save.txt',\
                                              rate_variation = False)
        self.MG94_IGC.update_by_x(x[:-1])
        self.MG94_IGC._loglikelihood()
        
        self.x= x
        # x = np.concatenate((self.MG94_IGC.x, np.log(tract_p)))

        self.IGC_sitewise_lnL_file = IGC_sitewise_lnL_file
        self.NOIGC_sitewise_lnL_file = NOIGC_sitewise_lnL_file
        self.MG94_IGC.get_sitewise_loglikelihood_summary(IGC_sitewise_lnL_file, False)
        self.MG94_IGC.get_sitewise_loglikelihood_summary(NOIGC_sitewise_lnL_file, True)
        
        
        self.outgroup_branch = [edge for edge in self.MG94_IGC.edge_list if edge[0] == 'N0' and edge[1] != 'N1'][0]
        Total_blen = sum([self.MG94_IGC.edge_to_blen[edge] for edge in self.MG94_IGC.edge_list if edge != self.outgroup_branch])

        #print 'a'
        self.hmmtract = HMMTract(IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file, State_List, Total_blen, self.MG94_IGC.tau, seq_index_file, model)
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
            self.save_x()
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
        #self.MG94_IGC.get_mle()
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

    def get_Hessian(self, One_Dimension = True):
        if One_Dimension:
            f = nd.Derivative(partial(self.hmmtract.objective_1D, False), n = 1)
            ff = nd.Derivative(partial(self.hmmtract.objective_1D, False), n = 2)
            grad = -f(self.hmmtract.x[1:])
            hess = -ff(self.hmmtract.x[1:])
        else:
            f = nd.Derivative(self._loglikelihood, n = 1)
            ff = nd.Hessdiag(self._loglikelihood)
            grad = -f(self.x)
            hess = -ff(self.x)
        return grad, hess

        
    

if __name__ == '__main__':
    pair = ["EDN", "ECP"]
    paralog = pair
    model = 'HKY'
    force = None
    nsites = None
    clock = None

    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
    newicktree = '../test/input_tree.newick'
    Force = None
    output_ctrl = ''
    summary_mat = []
    


    x = np.array([-6.974760482698901809e-01,
                  -5.370097148062192849e-01,
                  -7.240186071588703420e-01,
                  7.278666576839283309e-01,
                  #-1.661538873324752974e-01,
                  -6.389150587062580877e-01,
                  -1.627854846878343364e+00,
                  -1.118866867842108759e+00,
                  -3.421887713565998190e+00,
                  -1.856026261560800084e+00,
                  -3.362580076646019211e+00,
                  -2.407034335517133528e+00,
                  -4.244048731094335558e+00,
                  -4.133421357939401020e+00,
                  #0.0#-2.6929935812559425#-2.5
                  0.0
#                  -3.0832851705618527
                  ])
    x_mle = np.array([-0.6974872 , -0.53696826, -0.72405219,  0.7278009 , -0.16637466,
       -0.63867808, -1.62767928, -1.11891194, -3.42235663, -1.85612631,
       -3.36291767, -2.407368  , -4.2448557 , -4.1330508 , -3.07197501])

    print 
    print '**' + '_'.join(paralog)+ '**', output_ctrl


    alignment_file = '../test/EDN_ECP_Cleaned.fasta'
    cdna = True
    rate_variation = True
    if rate_variation:
        save_file = '../test/save/HMMJS_' + '_'.join(paralog) + '_' + model + '_rv_nonclock_save.txt'
    else:
        save_file = '../test/save/HMMJS_' + '_'.join(paralog) + '_' + model + '_nonclock_save.txt'    
    seq_index_file = '../test/' + '_'.join(paralog) +'_seq_index.txt'
        
    summary_path = '../test/Summary/'
    IGC_sitewise_lnL_file = '../test/Summary/' + '_'.join(paralog) + '_' + model + '_nonclock_sw_lnL.txt'
    NOIGC_sitewise_lnL_file = '../test/Summary/NOIGC_' + '_'.join(paralog) + '_' + model + '_nonclock_sw_lnL.txt'
    save_path = '../test/save/'

    seq_index_file = '../test/' + '_'.join(paralog) + '_seq_index.txt'

    if rate_variation:
        Ind = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = None, save_name = '../test/save/' + '_'.join(paralog) +'_' + model + '_rv_nonclock_save.txt', rate_variation = rate_variation)
    else:
        Ind = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = None, save_name = '../test/save/' + '_'.join(paralog) +'_' + model + '_nonclock_save.txt', rate_variation = rate_variation)

    Ind.get_mle()

    x = np.concatenate((Ind.x, [0.0]))

    test = HMMJSGeneconv(save_file, newicktree, alignment_file, paralog, summary_path, x, save_path, IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file,
                         state_list, seq_index_file, model,
                         rate_variation = rate_variation)
    test.update_by_x(x)

    self = test

    summary_file_1D = '../test/Summary/HMM_' + '_'.join(paralog) + '_' + model + '_nonclock_1D_summary.txt'
    summary_file_all_Dimension = '../test/Summary/HMM_' + '_'.join(paralog) + '_' + model + '_nonclock_all_summary.txt'

    log_p_list = np.log(3.0/np.array(range(3, 501)))
    plot_file = '../test/plot/HMM_' + '_'.join(paralog) + '_' + model + '_lnL_1D_surface.txt'
    test.plot_tract_p(log_p_list, plot_file)
    test.get_mle(display = True, two_step = True, One_Dimension = True)
    test.get_summary(summary_file_1D)
    #test.get_mle(display = True, two_step = False, One_Dimension = False)
    #test.get_summary(summary_file_all_Dimension)

    hess = test.get_Hessian(True)
    

######################################################################################
######################################################################################
######################################################################################
##    
##    gene_to_orlg_file = '../test/YDR418W_YEL054C_GeneToOrlg.txt'
##    alignment_file = '../test/YDR418W_YEL054C_Simulation.fasta'
##    paralog=['YDR418W', 'YEL054C']
##    model = 'MG94'
##    force = None
##    nsites = None
##    clock = None
##    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
##    save_file = '../test/save/HMMJS_' + '_'.join(paralog) + '_MG94_nonclock_save.txt'
##    summary_path = '../test/Summary/'
##    save_path = '../test/save/'
##    IGC_sitewise_lnL_file = '../test/Summary/Simulation_test_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
##    NOIGC_sitewise_lnL_file = '../test/Summary/Simulation_test_NOIGC_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
##
##    x_mle = np.array([ -0.84490517,  -0.72519025,  -1.74714445,   2.62465942,
##        -0.42422119,   1.87579805,  -4.54256365,  -9.        ,
##        -4.29195431,  -4.23903611,  -6.02963425,  -6.08910566,
##        -5.32986981,  -5.01906984, -17.38161318,  -5.31747459,
##        -5.35151234,  -6.04604204, 0.0])
##
##
##    newicktree = '../test/YeastTree.newick'
##    DupLosList = '../test/YeastTestDupLost.txt'
##    Force = None
##    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
##    node_to_pos = {'D1':0}
##    seq_index_file = '../test/YDR418W_YEL054C_seq_index.txt'
##
####    Ind_MG94_IGC = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = '../test/save/',
####                             save_name = '../test/save/simulation_test_save.txt')
####    Ind_MG94_IGC.get_mle(True, True, 0, 'BFGS')
##
##    test = HMMJSGeneconv(save_file, newicktree, alignment_file, paralog, summary_path, x_mle, save_path, IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file,
##                         state_list, seq_index_file)
##    test.update_by_x(x_mle)
##    test.get_mle(display = True, two_step = True, One_Dimension = True)
##
##    log_p_list = np.log(3.0/np.array(range(3, 501)))
##    plot_file = '../test/plot/Simulation_test_HMM_' + '_'.join(paralog) + '_lnL_1D_surface.txt'
##    test.plot_tract_p(log_p_list, plot_file)

    

    


