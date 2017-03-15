# A separate file for Simulator
# This simulation should simulate multigene family evolution with
# point mutation and IGC(interlocus gene conversion) processes
# along a given species tree with gene duplication loss history and arbitrary family size
# This simulator only simulates cdna for now
# Xiang Ji
# xji3@ncsu.edu

from IGCTractModel import IGCTractModel
from Tree import Tree
from PMModel import PMModel
from Common import divide_configuration, draw_from_distribution
import numpy as np
import cPickle, os
from copy import deepcopy
from math import floor


class Simulator:

    def __init__(self, pm_model_name, x_pm, rate_variation,   # PM Model
                 x_IGC, pm_IGC,                               # IGC Tract Model
                 tree_newick, DupLosList, x_rates,            # Tree input
                 terminal_node_list, node_to_pos,             # Configuration input
                 gene_to_orlg_file, seq_file, log_file,       # output info
                 seed_file,                                   # random seed file
                 seq_index_file,                              # sequence index file
                 cdna = True                                  # simulate cdna sequence only or not
                 ):

        self.tree = Tree(tree_newick, DupLosList, terminal_node_list, node_to_pos)
        self.x_rates = x_rates
        self.node_to_pos = node_to_pos
        self.terminal_node_list = terminal_node_list                 
        self.root_by_dup = self.tree.is_duplication_node(self.tree.phylo_tree.root.name)
        self.PMModel = PMModel(pm_model_name, x_pm, rate_variation)
        self.IGCModel = IGCTractModel(x_IGC, range(self.tree.n_orlg), pm_IGC)
        self.seq_index_file = seq_index_file # seq_index file location
        self.seq_index = None                # store sequence index information
        self.nsites    = None                # number of nucleotide to simulate, must agree with the overall span seq_index
        self.cdna      = cdna

        self.seq_file  = seq_file            # output sequence file
        self.log_file  = log_file            # output log file
        self.seed_file = seed_file           # random seed file

        # {node:{orthologous_group:sequence}}
        self.node_to_seq = dict()            # dictionary to store sequence at each node
        

        self.initiate()

    def __str__(self):  # overide for print function
        print self.PMModel
        print self.IGCModel
        print self.tree

        return 'IGC simulator output seq: ' + self.seq_file + '\n' + \
               'log file: ' + self.log_file + '\n'

    def initiate(self):
        self.get_seed_file()
        self.unpack_x_rates(self.x_rates)
        self.read_seq_index_file()

    def unpack_x_rates(self, x_rate):  # copied from PSJSGeneconv.py
        if self.root_by_dup:
            self.tree.unpack_x_rates(x_rate)
        else: # add additional constrain on the brach from root to first duplication node
            # The first element in x_rate is the combined rate for both N0 related edges
            # First duplication node should be on one of the two N0 related edges
            # Constraint: N0_Outgroup = 9 * N0_D1
            total_rate = np.exp(x_rate[0])
            outgroup_rate = x_rate[0] + np.log(0.9)
            N0_D1_rate = x_rate[0] + np.log(0.1)
            translated_rate = []
            rate_iter = 1
            for edge in self.tree.edge_list:
                if 'N0' in edge:
                    if 'D1' in edge:
                        translated_rate.append(N0_D1_rate)
                    else:
                        translated_rate.append(outgroup_rate)
                else:
                    translated_rate.append(x_rate[rate_iter])
                    rate_iter += 1
            self.tree.unpack_x_rates(translated_rate)        

    def get_seed_file(self):
        if os.path.isfile(self.seed_file):
            prng = cPickle.load(open(self.seed_file, 'r'))
            np.random.set_state(prng.get_state())
        else:
            prng = np.random.RandomState()
            cPickle.dump(prng, open(self.seed_file, 'w+'))            

    def read_seq_index_file(self):
        # The index should have columns:
        # nt_index, codon #, codon site for coding sequence
        # nt_index,  -1,  -1 for noncoding sequence
        seq_index = np.loadtxt(self.seq_index_file, dtype = int)
        if not self.cdna:
            seq_index[:, 1:] = -1
        else:
            # if it is protein coding sequence, then the length must be divisible by 3
            assert(seq_index[-1][0] %3 == 0 and len(seq_index) %3 == 0) 
            
        self.seq_index = seq_index
        assert(self.seq_index[0][0] == 1) # position start from 1
        self.nsites = self.seq_index[-1][0] # simulate the total overspan of the sequence

    def sim_root(self):
        root_name = self.tree.phylo_tree.root.name
        root_conf = self.tree.node_to_conf[root_name]
        root_orlg = divide_configuration(root_conf)
        self.node_to_seq[root_name] = dict()

        if self.PMModel.name == 'HKY':
            distn = [self.PMModel.parameters['Pi_' + nt] for nt in 'ACGT']
            for orlg in root_orlg['loc']:
                seq = draw_from_distribution(distn, self.nsites, 'ACGT')
                self.node_to_seq[root_name][orlg] = seq

    def sim(self):
        self.sim_root()
        # simulate sequences on each node by the path to each node in terminal_node_list


    def simulate_for_path(self, terminal_node_name):
        terminal_clade = self.tree.find_clade(terminal_node_name)
        path = self.tree.phylo_tree.get_path(terminal_clade)
        for clade in path:
            if clade.name in self.node_to_seq:
                continue
            else:
                father_clade = self.tree.find_parent_clade(clade.name)
                if father_clade.name in self.node_to_seq:
                    edge = (father_clade.name, terminal_node_name)
                    blen = self.tree.edge_to_blen[edge]
                    print
                else:
                    print 'The node cannot be recognised!'
                    assert(False)
            
    def sim_one_branch(self, starting_seq, blen, conf):
        branch_orlg = divide_configuration(conf)
        
        current_seq = deepcopy(starting_seq)
        branch_orlg = divide_configuration(conf)
        assert(all([orlg in branch_orlg['loc'] for orlg in starting_seq.keys()]))

        # Get sub IGC init matrix from
        ordered_orlg = sorted(branch_orlg['loc'])
        if len(ordered_orlg) == 1:
            Total_IGC_init_rate = 0.0
        else:
            branch_IGC_init_Q = np.zeros((len(ordered_orlg), len(ordered_orlg)), dtype = np.floating)
            for i in range(len(ordered_orlg)):
                for j in range(len(ordered_orlg)):
                    branch_IGC_init_Q[i, j] = self.IGCModel.Q_init[ordered_orlg[i], ordered_orlg[j]]
            
            IGC_init_rate_diag = branch_IGC_init_Q.sum(axis = 1) # row sum
            Total_IGC_init_rate = sum(IGC_init_rate_diag) * (self.nsites - 1 + 1/2)

        cummulate_time = 0.0

        while(cummulate_time < blen):
            # Now sample exponential distributed waiting time for next event
            # point mutation or IGC event
            # need to update total point mutation rate with every new event
            # no need to update total IGC rate on the same branch since it's modeled as context independent

            seq_rate_dict, Total_PM_rate = self.get_mutation_rate( current_seq )
            Total_rate = Total_PM_rate + Total_IGC_init_rate

            cummulate_time += np.random.exponential(1.0 / Total_rate)
            print cummulate_time

            if cummulate_time > blen :
                break
            else:
                # Now decide whether it's a point mutation or IGC event
                event = draw_from_distribution(np.array([Total_PM_rate, Total_IGC_init_rate]) / Total_rate,
                                               1, range(2))

                if event == 0:
                    # It's a point mutation event
                    self.get_one_point_mutation(current_seq, seq_rate_dict)
                elif event == 1:
                    # It's an IGC event
                    print

    def get_mutation_rate(self, seq): # modified from IGCSimulation.IGCSimulator
        #print self.current_seq
        # TODO: add in rv
        poisson_rate_sum = 0.0
        PM_diag_rates = self.PMModel.Q_mut.sum(axis = 1)

        seq_rate_dict = dict()
        for orlg in seq.keys():
            seq_rate = [PM_diag_rates['ACGT'.index(seq[orlg][i])] for i in range(self.nsites)]
            seq_rate_dict[orlg] = seq_rate
            poisson_rate_sum += sum(seq_rate)

        return seq_rate_dict, poisson_rate_sum


    def get_one_point_mutation(self, seq, seq_rate_dict): # modified from IGCSimulation.IGCSimulator
        # only allow all sequences to be of same length
        assert(len(set([len(seq_rate_dict[orlg]) for orlg in seq_rate_dict])) == 1)
        orlg_group = sorted(seq_rate_dict.keys())
        # Now concatenate all rates to one giant list to draw from distribution
        # Use a dumb but safe way rather than comprehension
        concatenated_rate = list()
        for orlg in orlg_group:
            concatenated_rate.extend(seq_rate_dict[orlg])
        concatenated_rate = np.array(concatenated_rate) / sum(concatenated_rate)

        # Now sample a point mutation position
        mut_pos = draw_from_distribution(concatenated_rate, 1, range(len(concatenated_rate)))
        mut_paralog_num = int(floor(mut_pos / len(seq_rate_dict[orlg_group[0]])))
        mut_paralog = orlg_group[mut_paralog_num]
        seq_pos = mut_pos - mut_paralog_num * len(seq_rate_dict[orlg_group[0]])

        # Now perform point mutation at the position
        old_state = seq[mut_paralog][seq_pos]
        prob = np.array(self.PMModel.Q_mut['ACGT'.index(old_state), :])
        new_state = 'ACGT'[draw_from_distribution(prob / sum(prob), 1, range(len(prob)))]
        seq[mut_paralog][seq_pos] = new_state

        # TODO: implement log
        # mutation_orlg, mut_pos, old_state, new_state
        mutation_info = [str(mut_paralog), str(seq_pos), old_state, new_state]
        print ' '.join(mutation_info)

        return seq, mutation_info

    def get_one_IGC_event(self, seq, sub_IGC_init_Q, sub_IGC_tract_Q, sub_total_IGC_init_Q, ordered_orlg):
        np.fill_diagonal(sub_IGC_tract_Q, 0.0)  # fill diagonal entries with 0, just in case..
        
        # sample an IGC event orlg pair        
        IGC_pos = draw_from_distribution(sub_total_IGC_init_Q.ravel(), 1, range(sub_total_IGC_init_Q.size))
        orlg_from_num = int(floor(IGC_pos / sub_total_IGC_init_Q.ndim))
        orlg_to_num = IGC_pos - orlg_from_num * sub_total_IGC_init_Q.ndim
        orlg_from = ordered_orlg[orlg_from_num]
        orlg_to   = ordered_orlg[orlg_to_num]

        # now sample a starting pos and tract length
        
        tract_p = sub_IGC_tract_Q[orlg_from_num, orlg_to_num]
        init_prob_array = np.array([1.0 / tract_p] + [1.0] * (self.nsites - 1))
        start_pos = draw_from_distribution(init_prob_array / sum(init_prob_array), 1, range(len(init_prob_array)))
        tract_length = np.random.geometric(tract_p, 1)
        stop_pos = start_pos + tract_length + 1
        seq, IGC_info = self.IGC_copy(start_pos, stop_pos, orlg_from, orlg_to, seq)

        return seq, IGC_info


    def IGC_copy(self, start_pos, stop_pos, orlg_from, orlg_to, seq):
        # make sure the two orlgs are in seq
        assert(orlg_from in seq and orlg_to in seq)
        template_seq = seq[orlg_from][start_pos:stop_pos]
        overide_seq = seq[orlg_to][start_pos:stop_pos]

        num_diff = sum([template_seq[i] != overide_seq[i] for i in range(len(template_seq))])

        # TODO: implement log
        IGC_info = [str(orlg_from), str(orlg_to), str(start_pos), str(stop_pos), str(num_diff)]
        print ' '.join(IGC_info)

        # Now perform the IGC
        print ''.join(template_seq), ''.join(overide_seq)
        for i in range(start_pos, stop_pos):
            seq[orlg_to][i] = seq[orlg_from][i]

        print ''.join(template_seq), ''.join(overide_seq)
        return seq, IGC_info
        
        

        
        

        
        
        
if __name__ == '__main__':
    gene_to_orlg_file = '../test/YDR418W_YEL054C_GeneToOrlg.txt'
    seq_file = '../test/YDR418W_YEL054C_Simulation.fasta'
    log_file = '../test/YDR418W_YEL054C_Simulation.log'
    seed_file = '../test/YDR418W_YEL054C_Simulation_seed.log'

    tree_newick = '../test/YeastTree.newick'
    DupLosList = '../test/YeastTestDupLost.txt'
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    seq_index_file = '../test/YDR418W_YEL054C_seq_index.txt'
    #nsites = 489

    pm_model_name = 'HKY'
    x_pm = np.log([0.4, 0.5, 0.2, 9.2])
    rate_variation = False

    x_IGC = [2.0, 0.3]
    init_pm = 'One rate'
    tract_pm = 'One rate'
    pm_IGC = [init_pm, tract_pm]

    x_rates = [-4.170654939766711422e+00,
               -5.674236262981605883e+00,
               -4.140979602575983520e+00,
               -4.344239699023852097e+00,
               -6.496123290482403334e+00,
               -6.063647134296714647e+00,
               -6.043806966727234276e+00,
               -5.111657692573940537e+00,
               -6.404488905061815451e+00,
               -5.467996717925044159e+00,
               -5.460686727891754799e+00,
               -6.459940982759793116e+00]
    
    test = Simulator(pm_model_name, x_pm, rate_variation,
                     x_IGC, pm_IGC, tree_newick, DupLosList, x_rates,
                     terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, log_file, seed_file, seq_index_file)

    self = test
    #print test
    test.sim()

    edge = ('N0', 'D1')
    edge = ('N0', 'kluyveri')
    blen = self.tree.edge_to_blen[edge]
    conf = self.tree.node_to_conf['N0']
    #conf = self.tree.node_to_conf['N2']
    starting_seq = self.node_to_seq['N0']

    current_seq = deepcopy(starting_seq)
    branch_orlg = divide_configuration(conf)
    #assert(all([orlg in branch_orlg['loc'] for orlg in starting_seq.keys()]))

    # Get sub IGC init matrix from
    ordered_orlg = sorted(branch_orlg['loc'])
    if len(ordered_orlg) == 1:
        Total_IGC_init_rate = 0.0
    else:
        branch_IGC_init_Q = np.zeros((len(ordered_orlg), len(ordered_orlg)), dtype = np.floating)
        branch_IGC_tract_Q = np.zeros((len(ordered_orlg), len(ordered_orlg)), dtype = np.floating)
        Total_IGC_init_Q = np.zeros((len(ordered_orlg), len(ordered_orlg)), dtype = np.floating)
        for i in range(len(ordered_orlg)):
            for j in range(len(ordered_orlg)):
                branch_IGC_init_Q[i, j] = self.IGCModel.Q_init[ordered_orlg[i], ordered_orlg[j]]
                branch_IGC_tract_Q[i, j] = self.IGCModel.Q_tract[ordered_orlg[i], ordered_orlg[j]]
                if i != j:
                    if branch_IGC_tract_Q[i, j] != 0:
                        Total_IGC_init_Q[i, j] = branch_IGC_init_Q[i, j] * (self.nsites - 1 + 1.0 / branch_IGC_tract_Q[i, j])
                
        IGC_init_rate_diag = branch_IGC_init_Q.sum(axis = 1) # row sum
        Total_IGC_init_rate = Total_IGC_init_Q.sum()

        sub_IGC_tract_Q = branch_IGC_tract_Q
        sub_total_IGC_init_Q = Total_IGC_init_Q
        sub_IGC_init_Q = branch_IGC_init_Q


    cummulate_time = 0.0

    while(cummulate_time < blen):
        # Now sample exponential distributed waiting time for next event
        # point mutation or IGC event
        # need to update total point mutation rate with every new event
        # no need to update total IGC rate on the same branch since it's modeled as context independent

        seq_rate_dict, Total_PM_rate = self.get_mutation_rate( current_seq )
        Total_rate = Total_PM_rate + Total_IGC_init_rate

        cummulate_time += np.random.exponential(1.0 / Total_rate)
        print cummulate_time

        

        if cummulate_time > blen :
            break
        else:
            # Now decide whether it's a point mutation or IGC event
            event = draw_from_distribution(np.array([Total_PM_rate, Total_IGC_init_rate]) / Total_rate,
                                           1, range(2))

            if event == 0:
                # It's a point mutation event
                current_seq, mutation_info = self.get_one_point_mutation(current_seq, seq_rate_dict)
            elif event == 1:
                # It's an IGC event
                current_seq, IGC_info = self.get_one_IGC_event(current_seq, branch_IGC_init_Q, branch_IGC_tract_Q, Total_IGC_init_Q, ordered_orlg)
        

    
    

    
    
    
    #print test.node_to_seq

    
