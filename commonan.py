# A separate file for Ancestral State Reconstruction
# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Xiang Ji
# xji3@ncsu.edu
# Tanchumin Xu
# txu7@ncsu.edu
from __future__ import print_function
import jsonctmctree.ll, jsonctmctree.interface
from JSGeneconv import JSGeneconv
from CodonGeneconv import *
from Func import *
from copy import deepcopy
import os
from Common import *
import numpy as np


def get_maxpro(list, nodecom):
    sites = np.zeros(shape=(len(nodecom), len(list)))
    for site in range(len(list)):
        i=0
        for node in nodecom:
            sites[i][site] = np.argmax(list[site][:, node])
            i=i+1
    return (sites)


class AncestralState:

    def __init__(self,
                 geneconv  # JSGeneconv analysis for now
                 ):

        self.geneconv = geneconv

    def get_mle(self):
        self.geneconv.get_mle()

    def get_scene(self):
        return self.geneconv.get_scene()

    def get_ancestral_state_dist(self):
        scene = self.get_scene()
        requests = [
            {"property": "DNDNODE"}
        ]
        # ref:https://jsonctmctree.readthedocs.io/en/latest/examples/yang2014/section_4_4_2_1_marginal/main.html?highlight=ancestral
        if isinstance(scene, list):  # ref: _loglikelihood() function in JSGeneconv.py
            for separate_scene in scene:
                j_in = {
                    'scene': separate_scene,
                    "requests": requests

                }
                j_out = jsonctmctree.interface.process_json_in(j_in)

        else:
            j_in = {
                'scene': scene,
                'requests': requests
            }
            j_out = jsonctmctree.interface.process_json_in(j_in)
            # TODO: finish this

        return j_out  # finish this too


if __name__ == '__main__':
    #######################
    ###########EDN ECP
    #######################lnL_scene = geneconv.get_scene()
    paralog = ['EDN', 'ECP']
    Force = None
    alignment_file = '../test/EDN_ECP_Cleaned.fasta'
    newicktree = '../test/input_tree.newick'
    Force = None
    model = 'MG94'
    save_name = '../test/save/Ind_MG94_EDN_ECP_nonclock_save.txt'

    geneconv = ReCodonGeneconv(newicktree, alignment_file, paralog, Model=model, Force=Force, clock=None,
                               save_path='../test/save/', save_name=save_name)
    test = AncestralState(geneconv)
    self = test
    j_out = test.get_ancestral_state_dist()
    aa = 0
    # for i in range(len(j_out["responses"][0][0])):
    #     print(j_out["responses"][0][0][i])
    #     aa=array(j_out["responses"][0][0][i])+aa
    # aa=self.get_scene()
    # print(aa["observed_data"])
    re = self.get_scene()
    list_for_iid = re["observed_data"]["iid_observations"]
    list_commonan = []
    # for i in range(len(list_for_iid)):
    for i in range(10):
        re["observed_data"]["iid_observations"] = [list_for_iid[i]]

        requests = [
            {"property": "DNDNODE"}
        ]
        j_in = {
            'scene': re,
            'requests': requests
        }
        j_out = jsonctmctree.interface.process_json_in(j_in)
        j_out_matrix = np.array(j_out["responses"][0][0])
        list_commonan.append(j_out_matrix)
    #  print(aa["process_definitions"][0]["row_states"])
    #  print(aa["process_definitions"][0]["column_states"])
    #  print(aa["process_definitions"][0]["transition_rates"])
    print(len(list_commonan))
    print(get_maxpro(list_commonan, [0,1,3,5,7]))