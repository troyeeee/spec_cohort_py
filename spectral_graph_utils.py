from os import curdir
from numpy.core.fromnumeric import sort
from numpy.core.shape_base import _block_setup
from scipy.sparse.linalg import eigsh
import numpy as np
from block import Block
import json

def cal_fielder_vector(mat):
    D = np.diag(np.sum(mat, axis= 0))
    L = np.matrix(D - mat)
    vals, vecs = eigsh(L, k=2, which='SM')
    fiedler_vector = vecs[:,[1]]
    # print("cal fielder vector, ", fiedler_vector)
    return fiedler_vector.tolist()

def linkage_status_from_vec_value_comp(vec):
    i=0
    res = []
    while(i<len(vec)):
        if vec[i][0] > vec[i+1][0]:
            res.append(i)
        else:
            res.append(i+1)
        i+=2
    return res

def linkage_to_haplotype(index):
    bases = '10'*len(index)
    return ''.join([bases[i] for i in index])

def init_weight_mat_with_d(mat_d):
    return np.zeros(shape=(2*mat_d, 2*mat_d), dtype='f')

def cal_haplotype(mat):
    fiedler_vector = cal_fielder_vector(mat)
    trivial_flag = True
    link = linkage_status_from_vec_value_comp(fiedler_vector)
    return link, trivial_flag, fiedler_vector

def cal_haplotype_in_block(blk_id, adj_mat):
    try:
        hap_res, flag, fiedler_vector = cal_haplotype(adj_mat)
    except Exception as e:
        print("error, no convergence, block id = ", blk_id)
        fdv = cal_fielder_vector(adj_mat)
        print(len(fdv), fdv)
        return [], False, ''
    if len(hap_res) == 0:
        print("error, link not trivialable and not all pos/neg, id = ", blk_id)
        fdv = cal_fielder_vector(adj_mat)
        return [], False, ''
    if not flag:
        print("error, link not trivialable, id = ", blk_id)
        fdv = cal_fielder_vector(adj_mat)
        return [], False, ''
    hap_seq = linkage_to_haplotype(hap_res)
    return hap_res,flag,hap_seq

def load_weight_by_pair(sorted_genotype_distance_map, genotype_seq_map, block, adj_mat):
    hom_w = 1
    for cur_idx in range(0, block.het_size):
        for next_idx in range(cur_idx+1, block.het_size):
            hom_flag = False
            major_allele_count = 0
            minor_allele_count = 0
            cur_pos_mat_id = 2*(cur_idx)
            cur_neg_mat_id = 2*(cur_idx)+1
            next_pos_mat_id = 2*(next_idx)
            next_neg_mat_id = 2*(next_idx)+1
            for distance,sorted_ref_ids in sorted_genotype_distance_map:
                if distance > 4:
                    break
                for ref_id in sorted_ref_ids:
                    ref_gt_seq = genotype_seq_map[ref_id]
                    cur_ref_gt = int(ref_gt_seq[cur_idx])
                    next_ref_gt = int(ref_gt_seq[next_idx])
                    if cur_ref_gt == 1 or next_ref_gt == 1:
                        continue
                    hom_flag = True
                    if cur_ref_gt == next_ref_gt:
                        adj_mat[cur_pos_mat_id][next_pos_mat_id]+=hom_w
                        adj_mat[cur_neg_mat_id][next_neg_mat_id]+=hom_w
                        adj_mat[next_pos_mat_id][cur_pos_mat_id]+=hom_w
                        adj_mat[next_neg_mat_id][cur_neg_mat_id]+=hom_w
                    else:
                        adj_mat[cur_pos_mat_id][next_neg_mat_id]+=hom_w
                        adj_mat[cur_neg_mat_id][next_pos_mat_id]+=hom_w
                        adj_mat[next_pos_mat_id][cur_neg_mat_id]+=hom_w
                        adj_mat[next_neg_mat_id][cur_pos_mat_id]+=hom_w
            if not hom_flag:
                    adj_mat[cur_pos_mat_id][next_pos_mat_id]+=hom_w/10
                    adj_mat[cur_neg_mat_id][next_neg_mat_id]+=hom_w/10
                    adj_mat[next_pos_mat_id][cur_pos_mat_id]+=hom_w/10
                    adj_mat[next_neg_mat_id][cur_neg_mat_id]+=hom_w/10
                    adj_mat[cur_pos_mat_id][next_neg_mat_id]+=hom_w/10
                    adj_mat[cur_neg_mat_id][next_pos_mat_id]+=hom_w/10
                    adj_mat[next_pos_mat_id][cur_neg_mat_id]+=hom_w/10
                    adj_mat[next_neg_mat_id][cur_pos_mat_id]+=hom_w/10
    return adj_mat

def phase_with_spec_graph_in_pair(sorted_genotype_distance_map, genotype_seq_map, block):
    res_l = []
    adj_mat = init_weight_mat_with_d(block.het_size)
    adj_mat = load_weight_by_pair(sorted_genotype_distance_map, genotype_seq_map, block, adj_mat)
    print(adj_mat)
    hap_res,flag,hap_seq = cal_haplotype_in_block(block.start, adj_mat)
    return hap_seq