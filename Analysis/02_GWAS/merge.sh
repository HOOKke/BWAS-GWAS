# merge 22 chromosome
plink \
  --bfile /mnt/data2/home02/khu/UKB_GWAS/CCA_combine/mental_sMRI/GWAS_IDP/snp_filter_missnp/ukb22828_c1_4\
  --merge-list /mnt/data2/home02/khu/UKB_GWAS/CCA_combine/mental_sMRI/GWAS_IDP/snp_filter_missnp/merge_list.txt \
  --make-bed \
  --out /mnt/data2/home02/khu/UKB_GWAS/CCA_combine/mental_sMRI/GWAS_IDP/GWAS/merge