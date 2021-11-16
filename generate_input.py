from pysam import VariantFile
import csv
from load_data import load_geno

def parse_vcf_to_file(path, hap_file, geno_file):
    bcf_in = VariantFile(path)
    h_file = open(hap_file, 'w')
    g_file = open(geno_file, 'w')
    samples = list(bcf_in.header.samples)
    hap_dict = dict()
    geno_dict = dict()
    for s in samples:
        hap_dict[s]=([],[])
        geno_dict[s]=[]
    for rec in bcf_in.fetch():
        for s in samples:
            hap_dict[s][0].append(str(rec.samples[s]['GT'][0]))
            hap_dict[s][1].append(str(rec.samples[s]['GT'][1]))
            geno_dict[s].append(str(rec.samples[s]['GT'][0] + rec.samples[s]['GT'][1]))
    for v in hap_dict.values():
        h_file.write(''.join(v[0])+'\n')
        h_file.write(''.join(v[1])+'\n')
    for v in geno_dict.values():
        g_file.write(''.join(v)+'\n')     
    bcf_in.close()
    h_file.close()
    g_file.close()

def load_poses_and_samples(file_path):
    bcf_in = VariantFile(file_path)
    poses_file = open("poses.txt", "w")
    samples_file = open("samples.txt", "w")
    for rec in bcf_in.fetch():
        poses_file.write(str(rec.pos) + "\n")
        # print(rec.pos)
    for rec in bcf_in.fetch():
        for sample in rec.samples.keys():
            samples_file.write(sample + "\n")
        break

if __name__== "__main__":
    # parse_vcf_to_file("EUR_test.vcf.gz", 'hap_file.txt', 'geno_file.txt')
    parse_vcf_to_file("ALL.chr22.poplit_cohort.vcf.gz", 'hap_file.txt', 'geno_file.txt')
    load_poses_and_samples("ALL.chr22.poplit_cohort.vcf.gz")
    hap_file = load_geno('hap_file.txt')
    geno_file = load_geno('geno_file.txt')