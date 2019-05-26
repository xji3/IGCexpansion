
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

def get_interior_node(list):
    n=list["node_count"]
    node=np.arange(n)
    interior_node=set(node)-set(re["observed_data"]["nodes"])
    c = [i for i in interior_node]
    return(c)

def translate_into_seq(promax,len_node,dict,model,len_se):
    list = []
    if model=='MG94':

        for i in range(len_node):
            p0="paralog0:"
            p1="paralog1:"
            for j in range(len_se):
                p0=p0+dict[(promax[i][j])//61]
                p1=p1+dict[(promax[i][j])%61]
            list.append(p0)
            list.append(p1)
    else:
        for i in range(len_node):
            p0="paralog0:"
            p1="paralog1:"
            for j in range(len_se):
                p0 = po+dict[(promax[i][j] )// 4]
                p1 = p1+dict[(promax[i][j] )% 4]
            list.append(p0)
            list.append(p1)

    return(list)

class AncestralState:

    def __init__(self,
                 geneconv  # JSGeneconv analysis for now
                 ):

        self.geneconv                 = geneconv
        self.ancestral_state_response = None
        self.scene                    = None
        self.num_to_state             = None
        self.num_to_node              = None

        if isinstance(geneconv, JSGeneconv):
            raise RuntimeError('Not yet implemented!')


    def get_mle(self):
        self.geneconv.get_mle()

    def get_scene(self):
        if self.scene is None:
            self.get_mle()
            self.scene = self.geneconv.get_scene()
        return self.scene

    def get_dict_trans(self):
        return self.geneconv.get_dict_trans()

    def get_ancestral_state_response(self):
        scene = self.get_scene()
        requests = [
            {'property': "DNDNODE"}
        ]
        # ref:https://jsonctmctree.readthedocs.io/en/latest/examples/yang2014/section_4_4_2_1_marginal/main.html
        if isinstance(scene, list):  # ref: _loglikelihood() function in JSGeneconv.py
            raise RuntimeError('Not yet tested.')
            separate_j_out = []
            for separate_scene in scene:
                j_in = {
                    'scene': separate_scene,
                    "requests": requests

                }
                j_out = jsonctmctree.interface.process_json_in(j_in)
                separate_j_out.append(j_out)
            result = separate_j_out

        else:
            j_in = {
                'scene': scene,
                'requests': requests
            }
            j_out = jsonctmctree.interface.process_json_in(j_in)
            if j_out['status'] is 'feasible':
                result = j_out['responses'][0]
            else:
                raise RuntimeError('Failed at obtaining ancestral state distributions.')
        return result

    def get_site_ancestral_dist(self, site_num):
        # site_num starts from 0 where 0 corresponds to the first site
        if self.ancestral_state_response is None:
            self.ancestral_state_response = self.get_ancestral_state_response()

        site_packed_dist = self.ancestral_state_response[site_num]
        node_state_prob_dict = {}
        for node_num in range(len(self.get_num_to_node())):
            state_prob_dict = {}
            for state_num in range(len(self.get_num_to_state())):
                node_state_prob = site_packed_dist[state_num][node_num]
                state_prob_dict[self.num_to_state[state_num]] = node_state_prob
            node_state_prob_dict[self.num_to_node[node_num]] =state_prob_dict
        return node_state_prob_dict

    def get_num_to_state(self):
        if self.num_to_state is None:
            if self.geneconv.Model == 'HKY':
                states = 'ACGT'
            elif self.geneconv.Model == 'MG94':
                states = geneconv.codon_nonstop
            self.num_to_state = {num:state for num, state in enumerate(product(states, repeat = 2))}                
        return self.num_to_state

    def get_num_to_node(self):
        if self.num_to_node is None:
            self.num_to_node = self.geneconv.num_to_node
        return self.num_to_node


if __name__ == '__main__':
    #######################
    ###########EDN ECP
    #######################
    paralog = ['EDN', 'ECP']
    Force = None
    alignment_file = '../test/EDN_ECP_Cleaned.fasta'
    newicktree = '../test/input_tree.newick'
    Force = None
    model = 'HKY'
    save_name = '../test/save/Ind_' + model + '_EDN_ECP_nonclock_save.txt'

    geneconv = ReCodonGeneconv(newicktree, alignment_file, paralog, Model=model, Force=Force, clock=None,
                               save_path='../test/save/', save_name=save_name)
    test = AncestralState(geneconv)
    self = test
    scene = test.get_scene()

    site_num = 0
    node_state_prob_dict = test.get_site_ancestral_dist(site_num)

##    aa = 0
    # for i in range(len(j_out["responses"][0][0])):
    #     print(j_out["responses"][0][0][i])
    #     aa=array(j_out["responses"][0][0][i])+aa
    # aa=self.get_scene()
    # print(aa["observed_data"])
##    re = self.get_scene()
##    list_for_iid = re["observed_data"]["iid_observations"]
##    list_commonan = []
##    for i in range(len(list_for_iid)):
##    # for i in range(3):
##        re["observed_data"]["iid_observations"] = [list_for_iid[i]]
##
##        requests = [
##            {"property": "DNDNODE"}
##        ]
##        j_in = {
##            'scene': re,
##            'requests': requests
##        }
##        j_out = jsonctmctree.interface.process_json_in(j_in)
##        j_out_matrix = np.array(j_out["responses"][0][0])
##        list_commonan.append(j_out_matrix)
##        # print(re["observed_data"]["iid_observations"])
##    #  print(aa["process_definitions"][0]["row_states"])
##    #  print(aa["process_definitions"][0]["column_states"])
##    #  print(aa["process_definitions"][0]["transition_rates"])
##    list_node=get_interior_node(re)
##    dict=self.get_dict_trans()
##    len_node=len(list_node)
##    len_se=len(list_commonan)
##    get_maxpro=get_maxpro(list_commonan,list_node)
##    # print(get_maxpro[2][2]%61)
##    translate=translate_into_seq(promax=get_maxpro,len_node=len_node,dict=dict,model=model,len_se=len_se)
##    print(translate)




