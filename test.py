from pysam import VariantFile

bcf_in = VariantFile("/home/yangshuo/Desktop/projects/data_repo/1000G/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", "r")

for rec in bcf_in.fetch():
    print(len(rec.samples.keys()))
    break