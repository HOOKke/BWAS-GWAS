for i in `seq 1 1 1`; do
  echo "$i"
  plink --bfile /mnt/data0/home/khu/UKB_GWAS/Gradient_neuroticism/sample_filter_mind/ukb_22828_c"$i"_2 \
  --keep /mnt/data2/home02/khu/UKB_GWAS/CCA_combine/mental_sMRI/GWAS_IDP/sample_filter_gene_mental/eid_gene_img.txt \
  --make-bed \
  --out /mnt/data2/home02/khu/UKB_GWAS/CCA_combine/mental_sMRI/GWAS_IDP/sample_filter_gene_mental/ukb22828_c"$i"_3
done