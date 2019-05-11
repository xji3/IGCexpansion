# A separate file for Ancestral State Reconstruction
# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Xiang Ji
# xji3@ncsu.edu
# Tanchumin Xu
# tchu7@ncsu.edu
from __future__ import print_function
import jsonctmctree.ll, jsonctmctree.interface
from JSGeneconv import JSGeneconv
from Func import *
from copy import deepcopy
import os
from Common import *

class AncestralState:

    def __init__(self,
                 geneconv                   # JSGeneconv analysis for now
                 ):

        self.geneconv = geneconv

    def get_mle(self):
        self.geneconv.get_mle()

    def get_scene(self):
        return self.geneconv.get_scene()

    def get_ancestral_state_dist(self):
        scene = self.get_scene()
        requests = {'property':'sndnode'}
        # ref:https://jsonctmctree.readthedocs.io/en/latest/examples/yang2014/section_4_4_2_1_marginal/main.html?highlight=ancestral
        if isinstance(scene, list):  # ref: _loglikelihood() function in JSGeneconv.py
            for separate_scene in scene:
                j_in = {
                    'scene' : separate_scene,
                    'requests' : request
                    }
                j_out = jsonctmctree.interface.process_json_in(j_in)
                # TODO: finish this
        else:
            j_in = {
                    'scene' : separate_scene,
                    'requests' : request
                    }
            j_out = jsonctmctree.interface.process_json_in(j_in)
            # TODO: finish this

        return None # finish this too
        
 


if __name__ == '__main__':
#######################
    ###########EDN ECP
#######################lnL_scene = geneconv.get_scene()
    gene_to_orlg_file = '../test/EDN_ECP_GeneToOrlg.txt'
    alignment_file = '../test/EDN_ECP_Cleaned_NewFormat.fasta'
    tree_newick = '../test/input_tree.newick'
    DupLosList = '../test/EDN_ECP_DupLost.txt'
    terminal_node_list = ['Chimpanzee', 'Gorilla', 'Orangutan', 'Macaque', 'Tamarin']
    node_to_pos = {'D1':0}
    seq_index_file = None
    pm_model = 'HKY'
    IGC_pm = 'One rate'
    space_list = None

    cdna = True
    allow_same_codon = True
    rate_variation = True
    save_file = '../test/save/JS_HKY_rv_EDN_ECP_nonclock_save.txt'
    summary_file = '../test/Summary/JS_HKY_rv_EDN_ECP_nonclock_summary.txt'
    x_js = np.log([ 0.4, 0.6, 0.7,  4.35588244, 0.8, 9.0,  0.3])
  
    force = None
    jsgeneconv = JSGeneconv(alignment_file, gene_to_orlg_file, cdna, tree_newick, DupLosList,x_js, pm_model, IGC_pm, rate_variation,
                      node_to_pos, terminal_node_list, save_file, force)

    test = AncestralState(jsgeneconv)
    self = test

    



