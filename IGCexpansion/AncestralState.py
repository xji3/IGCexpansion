# A separate file for Ancestral State Reconstruction
# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Xiang Ji
# xji3@ncsu.edu
# Tanchumin Xu
# tchu7@ncsu.edu
from __future__ import print_function
import jsonctmctree.ll, jsonctmctree.interface
from JSGeneconv import JSGeneconv
from CodonGeneconv import ReCodonGeneconv
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
        self.geneconv.get_mle() # make sure it's optimized first
        scene = self.get_scene()
        requests = [{'property':'DNDNODE'}]
        # ref:https://jsonctmctree.readthedocs.io/en/latest/examples/yang2014/section_4_4_2_1_marginal/main.html
        if isinstance(scene, list):  # ref: _loglikelihood() function in JSGeneconv.py
            for separate_scene in scene:
                j_in = {
                    'scene' : separate_scene,
                    'requests' : requests
                    }
                j_out = jsonctmctree.interface.process_json_in(j_in)
                # TODO: finish this
        else:
            j_in = {
                    'scene' : scene,
                    'requests' : requests
                    }
            j_out = jsonctmctree.interface.process_json_in(j_in)
            # TODO: finish this

        return j_out # finish this too
        
 


if __name__ == '__main__':
#######################
    ###########EDN ECP
#######################lnL_scene = geneconv.get_scene()
    paralog = ['EDN', 'ECP']
    Force = None
    alignment_file = '../test/EDN_ECP_Cleaned.fasta'
    newicktree = '../test/input_tree.newick'
    Force = None
    model = 'HKY'
    save_name = '../test/save/Ind_HKY_EDN_ECP_nonclock_save.txt'

    geneconv = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = None, save_path = '../test/save/', save_name = save_name)
    test = AncestralState(geneconv)
    self = test
    j_out = test.get_ancestral_state_dist()


    



