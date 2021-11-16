from block import Block
from pysam import VariantFile
import os

def load_list_from_file(list_file_path, int_flag):
    res_l = []
    with open(list_file_path,"r") as in_f:
        for i in in_f:
            if int_flag:
                res_l.append(int(i.strip()))
            else:
                res_l.append(i.strip())
    return res_l

def gen_final_hap(proband_blk, poses_file_path):
    final_hap_map = {}
    with open(poses_file_path, "r") as in_f:
        for idx,pos in enumerate(in_f):
            if idx in proband_blk.het_idxes:
                pos = int(pos.strip())
                het_idx = proband_blk.het_idxes.index(idx)
                het_gt = (int(proband_blk.het_pos_seq[het_idx]),int(proband_blk.het_neg_seq[het_idx]))
                final_hap_map[pos] = het_gt
    return final_hap_map

def regen_vcf_file(final_hap_map, temp_vcf_path, out_vcf_path):
    bcf_in = VariantFile(temp_vcf_path)
    bcf_out = VariantFile(out_vcf_path, 'w', header=bcf_in.header)
    i=0
    for idx,rec in enumerate(bcf_in.fetch()):
        if int(rec.pos) in final_hap_map.keys():
            if final_hap_map[int(rec.pos)] == (6,6):
                rec.samples[0]['GT'] = (1,0)
                rec.samples[0].phased = True
                print(idx, rec.samples[0]['GT'])
                bcf_out.write(rec)
                continue
            rec.samples[0]['GT'] = final_hap_map[int(rec.pos)]
            # print(rec.samples[0]['GT'], final_hap_map[int(rec.pos)])
            rec.samples[0].phased = True
            i=i+1
        if (rec.samples[0]['GT'][0]==rec.samples[0]['GT'][1]):
            rec.samples[0].phased = True
        bcf_out.write(rec)
    print(i,len(final_hap_map))


def gen_phased_vcf_pipeline(proband_blk,proband_idx,poses_file_path,sample_list_path):
    project_base_dir = "/home/yangshuo/Desktop/poplit/spec_cohort_py"
    cohort_file_path = project_base_dir + "/ALL.chr22.poplit_cohort.vcf.gz"
    output_file_dir_pattern = project_base_dir + "/output/phased_{}"
    sample_list = load_list_from_file(sample_list_path, False)
    final_hap_map = gen_final_hap(proband_blk,poses_file_path)
    sample_name = sample_list[proband_idx]
    output_file_dir = output_file_dir_pattern.format(sample_name)
    if not os.path.isdir(output_file_dir):
        os.makedirs(output_file_dir)
    tmp_vcf_path = output_file_dir + "/ALL.chr22.unphased.{}.vcf.gz".format(sample_name)
    unzip_tmp_vcf_path = output_file_dir + "/ALL.chr22.unphased.{}.vcf".format(sample_name)
    out_vcf_path = output_file_dir + "/ALL.chr22.phased.{}.vcf.gz".format(sample_name)
    unzip_out_vcf_path = output_file_dir + "/ALL.chr22.phased.{}.vcf".format(sample_name)
    stat_file_path = output_file_dir + "/stat_poplit.txt"
    # os cmd patterns
    extract_sample_cmd = "bcftools view -s {} -O z -o {} {}".format(sample_name, tmp_vcf_path, cohort_file_path)
    tabix_cmd = "tabix {}".format(tmp_vcf_path)
    gunzip_cmd_pattern = "gunzip {}"
    cal_err_cmd = "python calculate_haplotype_statistics.py -v1 {} -v2 {} > {}".format(unzip_out_vcf_path, unzip_tmp_vcf_path, stat_file_path)
    print(extract_sample_cmd)
    os.system(extract_sample_cmd)
    print(tabix_cmd)
    os.system(tabix_cmd)
    regen_vcf_file(final_hap_map, tmp_vcf_path, out_vcf_path)
    print(gunzip_cmd_pattern.format(tmp_vcf_path))
    os.system(gunzip_cmd_pattern.format(tmp_vcf_path))
    print(gunzip_cmd_pattern.format(out_vcf_path))
    os.system(gunzip_cmd_pattern.format(out_vcf_path))
    print(cal_err_cmd)
    os.system(cal_err_cmd)