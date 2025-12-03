library("TwoSampleMR")
# library("ieugwasr")
# library("plinkbinr")
# get_plink_exe()



# 2.读取结局GWAS数据，提取对应工具变量
out_dat <- read_outcome_data(
  # snps = exp_dat$SNP,
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


# 1.读取暴露GWAS数据，选取工具变量
# 1.1 读取数据
# exp_dat <- extract_instruments(outcomes='ieu-b-5102', access_token = NULL, p1=5*10^-8, clump=TRUE, r2=0.01)
# ieu-b-5102,ieu-b-42,
while (!exists("exp_dat")) {
  try({
    ## 提取暴露变量
    exp_dat <- extract_instruments(outcomes='ieu-b-102', access_token = NULL, p1=5*10^-8, clump=TRUE, r2=0.01)
  })
  Sys.sleep(2)
}
filtered_out_dat <- out_dat[out_dat$SNP %in% exp_dat$SNP, ]

gc()


# 3.对数据进行预处理，使其效应等位与效应量保持统一，这一步是必须的
dat <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = filtered_out_dat
)
rm(exp_dat)

# 4.MR分析与可视化
# 4.1 MR分析
res <- mr(dat, method_list=c(
  "mr_two_sample_ml",
  "mr_simple_median",
  "mr_weighted_median",
  "mr_ivw_mre",
  "mr_ivw_fe",
  "mr_ivw_radial"
))
res
res <- mr(dat, method_list=c(
  "mr_two_sample_ml",
  "mr_simple_median",
  "mr_weighted_median",
  "mr_ivw_mre",
  "mr_ivw_fe"
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
