from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
from IGCexpansion.JSGeneconv import JSGeneconv
import argparse, os
import numpy as np

def read_terminal_nodes(file_name):
    assert(os.path.isfile(file_name))  # make sure the file exists first
    with open(file, 'r') as f:
        terminal_node_list = [tip.strip() for tip in f.read().splitlines()]
    return terminal_node_list


if __name__ == '__main__':
    paralog1 = 'EDN'
    paralog2 = 'ECP'
    paralog = [paralog1, paralog2]
    rate_variation = True

    gene_to_orlg_file = './' + '_'.join(paralog) + '_GeneToOrlg.txt'
    alignment_file =  './' + '_'.join(paralog) + '_Cleaned_NewFormat.fasta'
    newicktree = './' + '_'.join(paralog) + '_tree.newick'
    DupLosList = './' + '_'.join(paralog) + '_DupLost.txt'
    Force = None
    terminal_node_list = ['Tamarin', 'Macaque', 'Orangutan', 'Gorilla', 'Chimpanzee']
    # terminal_node_list = read_terminal_nodes(args.terminal_node_file)
    node_to_pos = {'D1':0}
    seq_index_file = './' + '_'.join(paralog) + '_seq_index.txt'


###### Now get HKY+PSJS-IGC estimates
    IGC_pm = 'One rate'
    pm_model = 'HKY'
    log_paralog_folder = './log/' + '_'.join(paralog)
    save_paralog_folder = './save/' + '_'.join(paralog)
    summary_paralog_folder = './summary/' + '_'.join(paralog)

    if not os.path.isdir(log_paralog_folder):
        os.mkdir(log_paralog_folder)
    if not os.path.isdir(save_paralog_folder):
        os.mkdir(save_paralog_folder)
    if not os.path.isdir(summary_paralog_folder):
        os.mkdir(summary_paralog_folder)

    if rate_variation:
        save_file = save_paralog_folder + '/JS_' + pm_model + '_' + '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_save.txt'
        log_file = log_paralog_folder + '/JS_' + pm_model + '_' + '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_log.txt'
        summary_file = summary_paralog_folder + '/JS_' + pm_model + '_' + '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])  # initial guess paralog value
    else:
        save_file = save_sim_folder + '/JS_' + pm_model + '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
        log_file = log_sim_folder + '/JS_' + pm_model + '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_log.txt'
        summary_file = summary_sim_folder + '/JS_' + pm_model + '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])  # initial guess paralog value

    

    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, True, newicktree, DupLosList, x_js, pm_model, IGC_pm,
                         rate_variation, node_to_pos, terminal_node_list, save_file)

    test_JS.get_mle()
    test_JS.get_individual_summary(summary_file)

