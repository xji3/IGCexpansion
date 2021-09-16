from CodonGeneconv import *

class JointAnalysis:
    auto_save_step = 2

    def __init__(self,
                 alignment_file_list,
                 tree_newick,
                 paralog_list,
                 Model = 'MG94',
                 IGC_Omega = None,
                 nnsites = None,
                 Force = None,
                 Shared = None,
                 save_path = './save/',
                 save_name = None,
                 post_dup = 'N1'):
        # first check if the length of all input lists agree
        assert (len(set([len(alignment_file_list), len(paralog_list)])) == 1)
        # doesn't allow IGC-specific omega for HKY model
        assert(Model == 'MG94' or IGC_Omega is None)
        self.save_path     = save_path
        self.Model         = Model
        self.IGC_Omega     = IGC_Omega
        self.paralog_list  = paralog_list
        self.x             = None
        if Shared is None:
            self.shared_parameters = []
        else:
            self.shared_parameters = Shared
        grand_save_name, individual_save_names = self.get_save_file_names(save_name)
        self.geneconv_list = [ReCodonGeneconv(tree_newick, alignment_file_list[i], paralog_list[i],
                                              Model, IGC_Omega, nnsites, False, Force, save_path, individual_save_names[i], post_dup)
                              for i in range(len(alignment_file_list))]
        self.save_name     = grand_save_name

        self.auto_save = 0
        self.initialize_x()

    def initialize_x(self):
        if os.path.isfile(self.save_name):
            self.initialize_by_save(self.save_name)
            print('Successfully loaded parameter value from ' + self.save_name)
        else:
            single_x = self.geneconv_list[0].x
            shared_x = [single_x[i] for i in self.shared_parameters]
            unique_x = [single_x[i] for i in range(len(single_x)) if not i in self.shared_parameters] * len(self.geneconv_list)
            self.x = np.array(unique_x + shared_x)
        self.update_by_x(self.x)

    def get_save_file_names(self, save_name):
        if len(self.shared_parameters):
            model_string = self.Model + '_withSharing'
        else:
            model_string = self.Model

        if save_name is None:
            if self.IGC_Omega is None:
                general_save_name = self.save_path + 'Joint_' + model_string + '_' + str(len(self.paralog_list)) + '_pairs_grand_save.txt'
            else:
                general_save_name = self.save_path + 'Joint_' + model_string + '_twoOmega_' + str(len(self.paralog_list)) + '_pairs_grand_save.txt'
        else:
            general_save_name = save_name



        names = []
        for paralog in self.paralog_list:
            single_save_name = general_save_name.replace(str(len(self.paralog_list)) + '_pairs', '_'.join(paralog))
            names.append(single_save_name)

        return general_save_name, names

    def check_x_dim(self):
        assert(len(self.geneconv_list) > 1)
        if self.shared_parameters is None:
            shared_dim = 0
        else:
            shared_dim = len(self.shared_parameters)
        x_dim = sum([len(geneconv.x) for geneconv in self.geneconv_list]) - (len(self.geneconv_list) - 1) * shared_dim
        assert(len(self.x) == x_dim)

    def combine_x(self, uniq_x, shared_x):
        uniq_idx = 0
        shared_idx = 0
        x = []
        for i in range(len(self.geneconv_list[0].x)):
            if i in self.shared_parameters:
                x.append(shared_x[shared_idx])
                shared_idx += 1
            else:
                x.append(uniq_x[uniq_idx])
                uniq_idx += 1
        return x

    def update_by_x(self, x):
        self.check_x_dim()
        self.x = np.array(x)
        uniq_dim = len(self.geneconv_list[0].x) - len(self.shared_parameters)
        shared_x = self.x[len(self.geneconv_list) * uniq_dim:]
        for geneconv_idx in range(len(self.geneconv_list)):
            geneconv = self.geneconv_list[geneconv_idx]
            uniq_x = self.x[geneconv_idx * uniq_dim : (geneconv_idx + 1) * uniq_dim]
            geneconv.update_by_x(self.combine_x(uniq_x, shared_x))

    def get_original_bounds(self):
        bnds = [(None, -0.05)] * 3
        bnds.extend([(None, None)] * (len(self.geneconv_list[0].x) - 3))
        return bnds


    def combine_bounds(self):
        individual_bnds = self.get_original_bounds()
        combined_bounds = [individual_bnds[idx] for idx in range(len(individual_bnds)) if not idx in self.shared_parameters] * len(self.paralog_list) \
                          + [individual_bnds[idx] for idx in range(len(individual_bnds)) if idx in self.shared_parameters]
        return combined_bounds


    def objective_and_gradient(self, x):
        self.update_by_x(x)
        individual_results = [geneconv.objective_and_gradient(False, geneconv.x) for geneconv in self.geneconv_list]
        f = sum([result[0] for result in individual_results])
        uniq_derivatives = np.concatenate([[result[1][idx] for idx in range(len(result[1])) if not idx in self.shared_parameters] for result in individual_results])
        shared_derivatives = [[result[1][idx] for idx in range(len(result[1])) if idx in self.shared_parameters] for result in individual_results]
        g = np.concatenate((uniq_derivatives, np.sum(shared_derivatives, axis = 0)))
        self.auto_save += 1
        if self.auto_save == self.auto_save_step:
            self.save_x()
            self.auto_save = 0

        print('log likelihood = ', f)
        print('Current x array = ', self.x)
        print('Derivatives = ', g)
        return f, g


    def get_mle(self):
        self.update_by_x(self.x)

        guess_x = self.x

        result = scipy.optimize.minimize(self.objective_and_gradient, guess_x, jac=True, method='L-BFGS-B', bounds=self.combine_bounds())
        print (result)
        self.save_x()
        return result

    def save_x(self):
        save = self.x
        save_file = self.save_name
        np.savetxt(save_file, save.T)

    def initialize_by_save(self, save_file):
        self.x = np.loadtxt(open(save_file, 'r'))
        self.update_by_x(self.x)

if __name__ == '__main__':
    paralog = ['YLR406C', 'YDL075W']
    Force = None
    alignment_file = '../test/YLR406C_YDL075W_test_input.fasta'
    newicktree = '../test/YeastTree.newick'

    paralog_list = [paralog, paralog]
    IGC_Omega = None
    Shared = [5]
    alignment_file_list = [alignment_file, alignment_file]
    Model = 'MG94'

    joint_analysis = JointAnalysis(alignment_file_list,  newicktree, paralog_list, Shared = Shared,
                                   IGC_Omega = IGC_Omega, Model = Model, Force = Force,
                                   save_path = '../test/save/')
    # print(joint_analysis.objective_and_gradient(joint_analysis.x))
    joint_analysis.get_mle()



