# install.packages('devtools')
# devtools::install_github('MRCIEU/TwoSampleMR')

# setwd('D:/')
# devtools::install_local('./MRInstruments.zip')
# devtools::install_local('./MRPRESSO.zip')
# devtools::install_local('./TwoSampleMR.zip')
library("TwoSampleMR")

# library("ieugwasr")
# library("plinkbinr")
# get_plink_exe()

library(devtools)
install_github("phenoscanner/phenoscanner")
library(phenoscanner)


# 1.读取暴露GWAS数据，选取工具变量
# 1.1 读取数据
exp_dat <- read_exposure_data(
  filename = "D:/A/UK_Biobank/IDP/test4/results/GWAS/all_IDP/results/mode2/combine_result.P1.qassoc",
  sep = " ",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "MAF",
  pval_col = "P"
)
exp_dat0 <- exp_dat[exp_dat$eaf.exposure>=0.01,]

# 1.2 保证工具变量和暴露强相关
exp_dat1 <- exp_dat0[exp_dat0$pval.exposure<5*10^-8,]

# 1.3 clumping
exp_dat2 <- clump_data(exp_dat1, clump_r2 = 0.01, pop = "EUR")


# iv_snp <- exp_dat2$SNP
# res <- phenoscanner(snpquery="rs10840293")
# head(res$results)


# # 1.3 本地clumping
# exp_dat1 <- ld_clump(
#   dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, id=exp_dat$id.exposure), 
#   plink_bin = "D:/R-4.1.3/library/plinkbinr/bin/plink_Windows.exe", 
#   #欧洲人群参考基因组位置
#   bfile = "D:/A/LDSC/1kg.v3/EUR",
#   clump_kb = 10000,
#   clump_r2 = 0.01,
#   clump_p = 5e-08,
#   pop = "EUR"
# )
# exp_dat <- exp_dat[which(exp_dat$SNP %in% exp_dat1$rsid),]


# 2.读取结局GWAS数据，提取对应工具变量
# rm(exp_dat1)
# gc()
# out_dat <- extract_outcome_data(
#   snps = exp_dat2$SNP,
#   outcomes = 'ebi-a-GCST90027158'
# )

while (!exists("out_dat")) {
  try({
    ## 提取暴露变量
    out_dat <- extract_outcome_data(
      snps = exp_dat2$SNP,
      outcomes = 'ieu-a-1000'
    )
  })
  Sys.sleep(2)
}
# MH:ieu-b-5102,ieu-b-42,ieu-b-4763,ieu-b-4833,finn-b-F5_GAD,ieu-a-1183,ukb-d-20447,ebi-a-GCST005904,ebi-a-GCST009979,ebi-a-GCST009984,ieu-a-802,ieu-a-1007
# Cog:ieu-a-1066,ieu-a-1068,ebi-a-GCST90029012,ukb-a-505,ukb-b-2709,ebi-a-GCST006572,ieu-b-4838,

# 3.对数据进行预处理，使其效应等位与效应量保持统一，这一步是必须的
dat <- harmonise_data(
  exposure_dat = exp_dat2,
  outcome_dat = out_dat
)
rm(out_dat)

# 4.MR分析与可视化
# 4.1 MR分析
# mr_method_list()
res <- mr(dat, method_list=c(
  "mr_two_sample_ml",
  "mr_simple_median",
  "mr_weighted_median",
  "mr_ivw_mre",
  "mr_ivw_fe",
  "mr_ivw_radial"
))
#  "mr_ivw_radial","mr_raps"
res
res <- mr(dat, method_list=c(
  "mr_two_sample_ml",
  "mr_simple_median",
  "mr_weighted_median",
  "mr_ivw_mre",
  "mr_ivw_fe"
))
res
res <- mr(dat, method_list=c(
  "mr_two_sample_ml",
  "mr_egger_regression_bootstrap",
  "mr_simple_median",
  "mr_weighted_median",
  "mr_penalised_weighted_median",
  "mr_ivw",
  "mr_ivw_mre",
  "mr_ivw_fe",
  "mr_weighted_mode_nome",
  "mr_ivw_radial"
  ))
#  "mr_ivw_radial","mr_raps"
res
res <- mr(dat, method_list=c(
  "mr_two_sample_ml",
  "mr_egger_regression_bootstrap",
  "mr_simple_median",
  "mr_weighted_median",
  "mr_penalised_weighted_median",
  "mr_ivw",
  "mr_ivw_mre",
  "mr_ivw_fe",
  "mr_weighted_mode_nome"
))
res
res <- mr(dat, method_list=c("mr_ivw_radial"))
res
res <- mr(dat, method_list=c(
  "mr_two_sample_ml",
  "mr_egger_regression",
  "mr_egger_regression_bootstrap",
  "mr_simple_median",
  "mr_weighted_median",
  "mr_penalised_weighted_median",
  "mr_ivw",
  "mr_ivw_mre",
  "mr_ivw_fe",
  "mr_simple_mode",
  "mr_weighted_mode",
  "mr_weighted_mode_nome",
  "mr_simple_mode_nome",
  "mr_ivw_radial",
  "mr_uwr"
))
res
 

# 4.2 敏感性分析
# 4.2.1 异质性检验
mr_heterogeneity(dat)
# 4.2.2 水平多效性检验
mr_pleiotropy_test(dat)
# 4.2.3 留一检验：逐步剔除SNP后观察剩余的稳定性，理想的是剔除后变化不大，这和我们的meta分析剔除法很相似。
res_loo <- mr_leaveoneout(dat)
mr_leaveoneout_plot(res_loo)

# 4.3 可视化
# 4.3.1 散点图
p1 <- mr_scatter_plot(res, dat)
p1
# 4.3.2 森林图
res_single <- mr_singlesnp(dat) # 单个SNP的结果
p2 <- mr_forest_plot(res_single)
p2
# 4.3.3 漏斗图
p3 <- mr_funnel_plot(res_single)
p3
