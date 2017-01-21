# A file to store common functions that avoid iteratively import itself
# Xiang Ji
# xji3@ncsu.edu
import itertools

def divide_configuration(configuration):
    ortho_group_to_pos = dict(extent = {}, distinct = [], loc = [])
    # extent positions that represent same paralog (by the same ortholog group number) have to be in the same state
    # distinct positions don't change states, thus only track positions
    for pos in range(len(configuration)):
        if configuration[pos][1] == 1: # extent
            ortho_group = configuration[pos][0]
            if ortho_group in ortho_group_to_pos['extent']:
                ortho_group_to_pos['extent'][ortho_group].append(pos)
            else:
                ortho_group_to_pos['extent'][ortho_group] = [pos]
                ortho_group_to_pos['loc'].append(ortho_group)
        elif configuration[pos][1] == 0: # distinct
            ortho_group_to_pos['distinct'].append(pos)

    return ortho_group_to_pos 

def get_accessible_orlg_pair(conf_list):
    accessible_orlg_pair = list() # only store pair in an upper triangular fashion
    for conf in conf_list:
        ortho_group_to_pos = divide_configuration(conf)
        accessible_orlg_pair.extend(itertools.combinations(sorted(ortho_group_to_pos['loc']), 2))

    return sorted(list(set(accessible_orlg_pair)))
