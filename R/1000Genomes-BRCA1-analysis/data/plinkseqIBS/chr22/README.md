Run the following commands to generate the `.ibs` file in this directory:

```
gsutil cp gs://genomics-public-data/1000-genomes/vcf/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf .
gzip ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf
./plinkseq-0.10/pseq ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz ibs-matrix > ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.ibs
```

