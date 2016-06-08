import numpy as np
import os

def gen_summary_file_name(summary_path, model, clock, force, Dir, Dis, gBGC):
    prefix_summary = summary_path + model
    if force:
        prefix_summary = prefix_summary + '_Force'

    if Dir:
        prefix_summary = prefix_summary + '_Dir'

    if gBGC:
        prefix_summary = prefix_summary + '_gBGC'

    if clock:
        suffix_summary = '_clock_summary.txt'
    else:
        suffix_summary = '_nonclock_summary.txt'

    summary_file = prefix_summary + '_oldest_ADH1C_Dis_' + Dis + suffix_summary
    return summary_file


def summary_from_ind(oldest_paralog_list, summary_path, model, clock, force, Dir, Dis, gBGC):
    summary_file = gen_summary_file_name(summary_path, model, clock, force, Dir, Dis, gBGC)
    summary_mat = []
    for oldest_paralog in oldest_paralog_list:
        ind_summary_file = summary_file.replace('ADH1C', oldest_paralog)
        if os.path.isfile(ind_summary_file):
            res = np.loadtxt(open(ind_summary_file, 'r'))
            if len(np.atleast_1d(res)) > 1:
                summary_mat.append(res.tolist())
                label = open(ind_summary_file, 'r').readlines()[-1][2:-1]

    summary_file = summary_file.replace('_oldest_ADH1C', '')
    t = np.matrix(summary_mat)
    header = ' '.join(oldest_paralog_list)
    footer = label
    np.savetxt(open(summary_file, 'w+'), t.T, delimiter = ' ', header = header, footer = footer)


if __name__ == '__main__':
    summary_path = './summary/'
    model = 'HKY'
    oldest_paralog_list = ['ADH1A', 'ADH1B', 'ADH1C']
    clock = False
    force = True
    Dir   = False
    gBGC  = False
    Dis   = 'Free'

    summary_from_ind(oldest_paralog_list, summary_path, model, clock, force, Dir, Dis, gBGC)
    
        

    
    
