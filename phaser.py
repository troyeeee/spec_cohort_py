import json
from load_data import load_geno
from block import Block
from distance_tools import sort_by_IBS_distance
from spectral_graph_utils import phase_with_spec_graph_in_pair
from common_utils import comp
from result_utils import gen_phased_vcf_pipeline

def init_block(start, end, geno):
    blk = Block()
    blk.start = start
    blk.end = end
    blk.whole_seq = geno
    blk.whole_size = end-start
    blk.het_seq = extract_het_from_snp(blk.whole_seq)[0]
    blk.het_idxes = extract_het_from_snp(blk.whole_seq)[1]
    blk.het_poses = [i+start for i in blk.het_idxes]
    blk.het_size = len(blk.het_poses)
    blk.het_pos_seq = '1'*blk.het_size
    blk.het_neg_seq = '0'*blk.het_size
    return blk

def extract_het_from_snp(snp_seq):
    res_idx = []
    res_seq = ''
    for idx, i in enumerate(snp_seq):
        if i == '1':
            res_idx.append(idx)
            res_seq+=i
    return res_seq, res_idx

def get_idxes_of_item(lst, item):
    return [i for i, x in enumerate(lst) if x == item]

def extract_geno_from_idx(single_rec, idxes):
    res= ''
    for idx, i in enumerate(single_rec):
        if idx in idxes:
            res+=i
    return res

def get_geno_refs_from_idxes(surrogate_ref, idxes):
    res = list(map(lambda x: extract_geno_from_idx(x, idxes), surrogate_ref))
    return res

def fetch_block_ref(start, end, ref_lst):
    return [x[start:end] for x in ref_lst]

def first_iter(geno_list):
    first_iter_haps = list()
    # load one block for test
    start = 0
    end = len(geno_list[0])
    # end = 32
    inblock_geno_refs = fetch_block_ref(start, end, geno_list)
    # phase sample in turn
    for idx, geno in enumerate(inblock_geno_refs):
        proband = init_block(start, end, geno)
        proband_idx = idx
        # get the rest geno
        surrogate_inblock_refs = list(filter(lambda x: inblock_geno_refs.index(x) \
             != idx, inblock_geno_refs))
        sur_inblk_prbd_het_refs = get_geno_refs_from_idxes(surrogate_inblock_refs, \
            proband.het_idxes)
        # print(sur_inblk_prbd_het_refs)
        # phasing here
        sorted_genotype_distance_map,genotype_seq_map,tar_het_idx_set = sort_by_IBS_distance(proband_idx,inblock_geno_refs)
        res_het_seq = phase_with_spec_graph_in_pair(sorted_genotype_distance_map, genotype_seq_map, proband)
        proband.het_pos_seq = res_het_seq
        proband.het_neg_seq = comp(res_het_seq)
        gen_phased_vcf_pipeline(proband, idx, "poses.txt", "samples.txt")
        # print(json.dumps(proband, default=lambda obj:obj.__dict__))
        # store phased proband to first_iter_haps
        first_iter_haps.append(proband)
        if idx > 500:
            break

if __name__=="__main__":
    geno_file = "geno_file.txt"
    geno_list = load_geno(geno_file)
    first_iter(geno_list)
    