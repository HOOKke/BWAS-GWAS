#!/bin/bash
for i in `seq 1 1 22`; do
  echo "$i"
  plink \
  --bfile /mnt/data2/home02/khu/UKB_GWAS/CCA_combine/mental_sMRI/GWAS_IDP/sample_filter_white_unrelated/ukb22828_c"$i"_3 \
  --exclude /mnt/data2/home02/khu/UKB_GWAS/CCA_combine/mental_sMRI/GWAS_IDP/merge/merge-merge.missnp \
  --make-bed \
  --out /mnt/data2/home02/khu/UKB_GWAS/CCA_combine/mental_sMRI/GWAS_IDP/snp_filter_missnp/ukb22828_c"$i"_4
done