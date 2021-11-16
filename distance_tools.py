from typing import OrderedDict

def sort_by_IBS_distance(tar_id, genotype_list):
    print(tar_id, len(genotype_list[tar_id]))
    genotype_distance_map = {}
    genotype_seq_map = {}
    tar_het_idx_set = set()
    tar_genotype = genotype_list[tar_id]
    for idx,genotype in enumerate(genotype_list):
        IBS_SCORE_DISTANCE = 0
        surrogate_gt_seq = ""
        if idx == tar_id:
            continue
        for snp_site in range(0,len(genotype_list[tar_id])):
            # print("snp site = ", snp_site)
            tar_snp_gt = int(tar_genotype[snp_site])
            ref_snp_gt = int(genotype[snp_site])
            if abs(ref_snp_gt-tar_snp_gt) == 2:
                IBS_SCORE_DISTANCE += 2
            if tar_snp_gt == 1:
                if idx == 0:
                    tar_het_idx_set.add(snp_site)
                surrogate_gt_seq += str(ref_snp_gt)
        genotype_seq_map[idx] = surrogate_gt_seq
        if IBS_SCORE_DISTANCE not in genotype_distance_map:
            genotype_distance_map[IBS_SCORE_DISTANCE] = []
        genotype_distance_map[IBS_SCORE_DISTANCE].append(idx)
    sorted_genotype_distance_map = sorted(genotype_distance_map.items(),key=lambda x:x[0])
    return sorted_genotype_distance_map,genotype_seq_map,tar_het_idx_set