# A separate file for joint estimation of HMM + JSIGC model
# Xiang Ji
# xji3@ncsu.edu

from HMMTract import *
from CodonGeneconv import ReCodonGeneconv

class HMMJSGeneconv:
    auto_save_step = 2
    def __init__(self, # First, ReCodonGeneconv objects
                 newicktree, alignment_file, paralog, summary_path, x,
                 # Now HMMTract objects
                 IGC_sitewise_lnL_file, Force_sitewise_lnL_file,
                 State_List, tau, seq_index_file,
                 model = 'MG94',
                 force = None, nsites = None, clock = None):
        
        self.MG94_IGC = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = force, clock = clock)
        self.MG94_Force = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = {5:0.0}, clock = clock)
        self.MG94_IGC.update_by_x(x)
        self.MG94_Force.update_by_x(x)
        if not os.path.isfile(IGC_sitewise_lnL_file):
            self.MG94_IGC.get_sitewise_loglikelihood_summary(IGC_sitewise_lnL_file)
        if not os.path.isfile(Force_sitewise_lnL_file):
            self.MG94_Force.get_sitewise_loglikelihood_summary(Force_sitewise_lnL_file)
        
        self.outgroup_branch = [edge for edge in self.MG94_IGC.edge_list if edge[0] == 'N0' and edge[1] != 'N1'][0]
        Total_blen = sum([self.MG94_IGC.edge_to_blen[edge] for edge in self.MG94_IGC.edge_list if edge != self.outgroup_branch])
        
        self.hmmtract = HMMTract(IGC_sitewise_lnL_file, Force_sitewise_lnL_file, State_List, Total_blen, tau, seq_index_file)


if __name__ == '__main__':
    pair = ["EDN", "ECP"]
    paralog = pair

    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
    newicktree = '../test/input_tree.newick'
    Force = None
    output_ctrl = ''
    summary_mat = []

    x = np.array([-6.972787845489714087e-01,
                  -5.371080077441253708e-01,
                  -7.240047380894862883e-01,
                  7.238578797129271436e-01,
                  -1.898373482480217311e-01,
                  -4.600735470664583104e-01,
                  -1.613942190237640295e+00,
                  -1.123019321880592836e+00,
                  -3.443903685022769334e+00,
                  -1.854709904478881510e+00,
                  -3.365339758069785692e+00,
                  -2.410226233318178313e+00,
                  -4.249147091385517605e+00,
                  -4.137233772266141862e+00])

    print 
    print '**' + '_'.join(paralog)+ '**', output_ctrl
    
    alignment_file = '../test/EDN_ECP_Cleaned.fasta'
    summary_path = '../test/Summary/'
    MG94_IGC = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
    
    IGC_sitewise_lnL_file = '../test/Summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
    Force_sitewise_lnL_file = '../test/Summary/Force_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
    Total_blen = sum([MG94_IGC.edge_to_blen[edge] for edge in MG94_IGC.edge_list if edge != ('N0', 'Tamarin')])
    seq_index_file = './' + '_'.join(paralog) + '_seq_index.txt'

    test = HMMJSGeneconv(newicktree, alignment_file, paralog, summary_path, x, IGC_sitewise_lnL_file, Force_sitewise_lnL_file,
                         state_list, 1.0, seq_index_file)
