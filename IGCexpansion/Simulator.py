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
                 gene_to_orlg_file, seq_file,                 # output info
                 IGC_log_file, PM_log_file,                   # log_files
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
        self.gene_to_orlg_file = gene_to_orlg_file          # input from data class which is not used here
        self.gene_to_orlg   = dict()         # dictionary used to store gene ortholog group info
        self.orlg_to_gene   = dict()         # dictionary used to store gene ortholog group info
        self.seq_index_file = seq_index_file # seq_index file location
        self.seq_index = None                # store sequence index information
        self.nsites    = None                # number of nucleotide to simulate, must agree with the overall span seq_index
        self.cdna      = cdna

        self.seq_file      = seq_file            # output sequence file
        self.IGC_log_file  = IGC_log_file        # IGC output log file
        self.PM_log_file   = PM_log_file         # PM output log file
        self.seed_file     = seed_file           # random seed file

        # {node:{orthologous_group:sequence}}
        self.node_to_seq = dict()            # dictionary to store sequence at each node
        

        self.initiate()

    def append_to_log_file(self, new_line, case):
        if case == 'IGC':
            with open(self.IGC_log_file, 'a') as g:
                g.write('\t'.join([str(item) for item in new_line]) + '\n')
        elif case == 'PM':
            with open(self.PM_log_file, 'a') as g:
                g.write('\t'.join([str(item) for item in new_line]) + '\n')
            

    def __str__(self):  # overide for print function
        print self.PMModel
        print self.IGCModel
        print self.tree

        return 'IGC simulator output seq: ' + self.seq_file + '\n' + \
               'IGC log file: ' + self.IGC_log_file + \
               'PM log file: ' + self.PM_log_file + '\n'

    

    def get_gene_to_orlg(self):  # copied from data class
        assert(os.path.isfile(self.gene_to_orlg_file))
        with open(self.gene_to_orlg_file, 'rb') as f:
            for line in f:
                items = line.split()
                if items:
                    gene = items[0]
                    orlg = int(items[1])
                    self.gene_to_orlg[gene] = orlg
                    if orlg in self.orlg_to_gene:
                        self.orlg_to_gene[orlg].append(gene.split('__'))
                    else:
                        self.orlg_to_gene[orlg] = [gene.split('__')]

    def initiate(self):
        self.get_seed_file()
        self.get_gene_to_orlg()
        self.unpack_x_rates(self.x_rates)
        self.read_seq_index_file()
        # Now update log files
        with open(self.IGC_log_file, 'w+') as f:
            f.write('\t'.join(['edge', 'time', 'orlg_from', 'orlg_to', 'start_pos', 'stop_pos', 'num_diff', 'template_seq', 'overide_seq']) + '\n')
        with open(self.PM_log_file, 'w+') as g:
            g.write('\t'.join(['edge', 'time', 'mut_orlg', 'mut_pos', 'old_state', 'new_state']) + '\n')
            

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

    def sim(self, display = False):
        self.sim_root()
        # simulate sequences on each node by the path to each node in terminal_node_list
        for terminal_node_name in self.terminal_node_list:
            self.simulate_for_path(terminal_node_name, display)


    def simulate_for_path(self, terminal_node_name, display):
        terminal_clade = self.tree.find_clade(terminal_node_name)
        path = self.tree.phylo_tree.get_path(terminal_clade)
        for clade in path:
            if clade.name in self.node_to_seq:
                continue
            else:
                father_clade = self.tree.find_parent_clade(clade.name)
                if father_clade.name in self.node_to_seq:
                    edge = (father_clade.name, clade.name)
                    self.sim_one_branch(edge, display)
                else:
                    print 'The node cannot be recognised!'
                    assert(False)
            
    def sim_one_branch(self, edge, display):
        # First, make sure this branch is not simulated
        assert(edge[0] in self.node_to_seq and not edge[1] in self.node_to_seq)
        blen = self.tree.edge_to_blen[edge]
        starting_seq = self.node_to_seq[edge[0]]
        conf = self.tree.node_to_conf[edge[0]]
        branch_orlg = divide_configuration(conf)
        
        current_seq = deepcopy(starting_seq)  # it's passed on to the next node, need a new allocation of memory
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
            if display:
                print cummulate_time
           

            if cummulate_time > blen :
                break
            else:
                # Now decide whether it's a point mutation or IGC event
                event = draw_from_distribution(np.array([Total_PM_rate, Total_IGC_init_rate]) / Total_rate,
                                               1, range(2))

                if event == 0:
                    # It's a point mutation event
                    current_seq, mutation_info = self.get_one_point_mutation(current_seq, seq_rate_dict, display)
                    to_write_info = ['_'.join(edge), str(cummulate_time)] + mutation_info
                    self.append_to_log_file(to_write_info, 'PM')
                    
                elif event == 1:
                    # It's an IGC event
                    current_seq, IGC_info = self.get_one_IGC_event(current_seq, branch_IGC_init_Q, branch_IGC_tract_Q, Total_IGC_init_Q, ordered_orlg, display)
                    to_write_info = ['_'.join(edge), str(cummulate_time)] + IGC_info
                    self.append_to_log_file(to_write_info, 'IGC')
                else:
                    # draw from distribution failure
                    assert(False)

        # Now need to pass the seq to new node
        # need to consider gene duplication loss events here
        self.pass_seq_to_node(current_seq, edge)

    def pass_seq_to_node(self, seq, edge):
        father_node = edge[0]
        child_node = edge[1]
        assert(self.tree.find_parent_clade(child_node).name == father_node) # make sure it's as expected

        father_orlg = divide_configuration(self.tree.node_to_conf[father_node])
        child_orlg = divide_configuration(self.tree.node_to_conf[child_node])

        if set(father_orlg['loc']) == set(child_orlg['loc']):
            self.node_to_seq[child_node] = seq
        else:
            father_changed_orlg = set(father_orlg['loc']) - set(child_orlg['loc'])
            child_changed_orlg = set(child_orlg['loc']) - set(father_orlg['loc'])
            assert(len(father_changed_orlg) == 1) # only allow one gene to duplicate or lose for now, for tandem duplication, should allow multiple
            if len(child_changed_orlg) == 2: # duplication event
                give_birth_orlg = list(father_changed_orlg)[0]
                assert(set(self.tree.dup_events[give_birth_orlg]) == set(child_changed_orlg))
                for new_orlg in child_changed_orlg:
                    seq[new_orlg] = deepcopy(seq[give_birth_orlg])
                seq.pop(give_birth_orlg)
            elif len(child_changed_orlg) == 0: # a gene loss event
                lost_orlg = list(father_changed_orlg)[0]
                seq.pop(lost_orlg)
            else:
                assert(False) # should not come to this case

            self.node_to_seq[child_node] = seq
 
    def get_mutation_rate(self, seq): # modified from IGCSimulation.IGCSimulator
        # print self.current_seq
        # TODO: add in rv
        poisson_rate_sum = 0.0
        PM_diag_rates = self.PMModel.Q_mut.sum(axis = 1)

        seq_rate_dict = dict()
        for orlg in seq.keys():
            seq_rate = [PM_diag_rates['ACGT'.index(seq[orlg][i])] for i in range(self.nsites)]
            # My lazy compromise
            if self.PMModel.rate_variation:
                for i in range(len(seq_rate) / 3):
                    seq_rate[3 * i + 1] = seq_rate[3 * i + 1] * self.PMModel.parameters['r2']
                    seq_rate[3 * i + 2] = seq_rate[3 * i + 2] * self.PMModel.parameters['r3']
            seq_rate_dict[orlg] = seq_rate
            poisson_rate_sum += sum(seq_rate)

        return seq_rate_dict, poisson_rate_sum


    def get_one_point_mutation(self, seq, seq_rate_dict, display): # modified from IGCSimulation.IGCSimulator
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
        if display:
            print ' '.join(mutation_info)

        return seq, mutation_info

    def get_one_IGC_event(self, seq, sub_IGC_init_Q, sub_IGC_tract_Q, sub_total_IGC_init_Q, ordered_orlg, display):
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
        tract_length = np.random.geometric(tract_p, 1)[0]
        stop_pos = start_pos + tract_length + 1
        if stop_pos > self.nsites:
            stop_pos = self.nsites
        #print start_pos, stop_pos, self.nsites
        seq, IGC_info = self.IGC_copy(start_pos, stop_pos, orlg_from, orlg_to, seq, display)

        return seq, IGC_info


    def IGC_copy(self, start_pos, stop_pos, orlg_from, orlg_to, seq, display):
        # make sure the two orlgs are in seq
        assert(orlg_from in seq and orlg_to in seq)
        template_seq = seq[orlg_from][start_pos:stop_pos]
        overide_seq = seq[orlg_to][start_pos:stop_pos]

        num_diff = sum([template_seq[i] != overide_seq[i] for i in range(len(template_seq))])

        # TODO: implement log
        IGC_info = [str(orlg_from), str(orlg_to), str(start_pos), str(stop_pos), str(num_diff)]
        if display:
            print ' '.join(IGC_info)

        # Now perform the IGC event
            print ''.join(template_seq), ''.join(overide_seq)

        IGC_info.extend([''.join(template_seq), ''.join(overide_seq)])
        for i in range(start_pos, stop_pos):
            seq[orlg_to][i] = seq[orlg_from][i]

        return seq, IGC_info


    def output_seq(self):
        # output sequence into seq_file using fasta format
        # use self.seq_index information to select sites being written into the file
        with open(self.seq_file, 'w+') as f:
            for name in self.terminal_node_list:
                for orlg in self.node_to_seq[name]:
                    gene_name_list = [i for i in self.orlg_to_gene[orlg] if i[0] == name ]
                    assert(len(gene_name_list) == 1)
                    f.write('>' + '__'.join(gene_name_list[0]) + '\n')
                    f.write(''.join(self.node_to_seq[name][orlg]) + '\n')
                    
            

        
        

        
        

        
        
        
if __name__ == '__main__':
    gene_to_orlg_file = '../test/YDR418W_YEL054C_GeneToOrlg.txt'
    seq_file = '../test/YDR418W_YEL054C_Simulation.fasta'
    IGC_log_file = '../test/YDR418W_YEL054C_Simulation_IGC.log'
    PM_log_file  = '../test/YDR418W_YEL054C_Simulation_PM.log'
    seed_file = '../test/YDR418W_YEL054C_Simulation_seed.log'

    tree_newick = '../test/YeastTree.newick'
    DupLosList = '../test/YeastTestDupLost.txt'
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    seq_index_file = '../test/YDR418W_YEL054C_seq_index.txt'
    #nsites = 489

    pm_model_name = 'HKY'
    x_pm = np.log([0.4, 0.5, 0.2, 9.2, 0.4, 5.0])
    rate_variation = True

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
                     terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, IGC_log_file, PM_log_file, seed_file, seq_index_file)

    self = test
    print test
    display = True
    test.sim(display = display)
    test.output_seq()

##    print "Still simulating?"
##    edge = ('N0', 'D1')
##    test.sim_one_branch(edge, True)
##    
##    edge = ('D1', 'N1')
##    test.sim_one_branch(edge, True)

    


    
    
    #print test.node_to_seq

    
