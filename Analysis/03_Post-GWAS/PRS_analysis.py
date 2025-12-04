import pandas as pd
import numpy as np
import scipy.stats as stats
import copy
import statsmodels.api as sm
import matplotlib.pyplot as plt
import os
import seaborn as sns
from scipy.stats.mstats import zscore


# # p_list = [0.5, 0.1, 0.05, 0.01, 1E-3, 1E-4, 1E-5, 1E-6, 1E-7, 5E-8]
# # for p in p_list:
# #     snp_list = './mode2/snp_list_' + str(p) + '.txt'
# #     print("plink "
# #           "--bfile /mnt/data2/home02/khu/UKB_GWAS/all_merge/GWAS/merge "
# #           "--score ./mode2/score.txt "
# #           "--extract " + snp_list + " --out ./mode2/PRS_" + str(p))
#
#
# # # prepare subject: IDP/mental/cog & gene, and not included in GWAS
# # GWAS_eid = pd.read_csv(r'D:\A\UK_Biobank\IDP\test4\results\GWAS\all_IDP\eid\eid_gwas_new_old_34869.csv')
# # gene_eid = pd.read_csv('gene_eid.csv')
# # rest_gene_eid = gene_eid[~gene_eid['eid'].isin(GWAS_eid['old'])]  # old 452426
# # eid_mapping = pd.read_csv(r'D:\A\UK_Biobank\ID_mapping\eid_mapping.csv')
# # overlap = eid_mapping[eid_mapping['old'].isin(rest_gene_eid['eid'])]  # new and old
# # # overlap.to_csv('rest_gene_eid.csv', index=False)
# # # IDP
# # IDP_eid = pd.read_csv(r'D:\A\UK_Biobank\IDP\test4\results\IDP\idp_sMRI.csv')  # 42301
# # IDP_gene = IDP_eid[IDP_eid['eid'].isin(overlap['new'])]
# # IDP_gene.to_csv('IDP_gene.csv', index=False)  # 6331
# # mental_gene = pd.read_csv('../correlation_mental/mental_gene.csv')  # 82209
# # IDP_gene = pd.read_csv('IDP_gene.csv')  # 6331
# # overlap = mental_gene[mental_gene['eid'].isin(IDP_gene['eid'])]  # 2656
# # mental
# # mental_eid = pd.read_csv('D:/A/UK_Biobank/IDP/test3/mental_combine_2/results/feature/mental_zscore.csv')  # 100550
# # mental_gene = mental_eid[mental_eid['eid'].isin(overlap['new'])]
# # mental_gene.to_csv('mental_gene.csv', index=False)  # 82209
# # cog
# # cog_eid = pd.read_csv('D:/A/UK_Biobank/IDP/test3/cognition/results/cognition/cognition_2_nonan_norm.csv')  # 35620
# # cog_gene = cog_eid[cog_eid['eid'].isin(overlap['new'])]
# # cog_gene.to_csv('./correlation_cog/cog_gene.csv', index=False)  # 12463
# # IQ
# # IQ_eid = pd.read_csv('D:/A/UK_Biobank/IDP/test3/cognition/results/cognition/cognition_2.csv')
# # IQ_eid = IQ_eid[['eid', '20016-2.0']]
# # IQ_eid = IQ_eid.dropna(axis=0, how='any')  # 51481
# # IQ_eid = IQ_eid[IQ_eid['eid'].isin(overlap['new'])]
# # IQ_eid.to_csv('IQ_gene.csv', index=False)
# # print(IQ_eid.shape[0])  # 17830
#
#
# # # prepare PRS: PRS & mental/cog
# # threshold = '0.5'
# # data = pd.read_csv('./correlation_cog/cog_gene.csv')
# # prs1 = pd.read_csv('./PRS/mode1/PRS_' + threshold + '.csv')
# # prs2 = pd.read_csv('./PRS/mode2/PRS_' + threshold + '.csv')
# # eid_mapping = pd.read_csv(r'D:\A\UK_Biobank\ID_mapping\eid_mapping.csv')
# # overlap1 = eid_mapping[eid_mapping['old'].isin(prs1['IID'])]
# # overlap2 = eid_mapping[eid_mapping['old'].isin(prs2['IID'])]
# # prs1.columns = ['old', 'score']
# # prs2.columns = ['old', 'score']
# # prs1 = pd.merge(overlap1, prs1, on='old')
# # prs2 = pd.merge(overlap2, prs2, on='old')
# # my_prs1 = prs1[prs1['new'].isin(data['eid'])]
# # my_prs2 = prs2[prs2['new'].isin(data['eid'])]
# # my_prs1 = my_prs1.sort_values(by='new')
# # my_prs2 = my_prs2.sort_values(by='new')
# # if not os.path.exists('./correlation_cog/PRS_' + threshold):
# #     os.makedirs('./correlation_cog/PRS_' + threshold)
# # my_prs1.to_csv('./correlation_cog/PRS_' + threshold + '/my_prs1.csv', index=False)
# # my_prs2.to_csv('./correlation_cog/PRS_' + threshold + '/my_prs2.csv', index=False)
#
#
# 删除mental health和prs离群值
def delete_outlier(my_data, col):
    s = my_data.iloc[:, col]
    threshold = s.std() * 3
    outliers = []
    for y in s:
        if np.abs(y - s.mean()) > threshold:
            outliers.append(y)
            index = my_data[my_data.iloc[:, col] == y].index.to_list()
            my_data.drop(index, inplace=True)


# correlation: PRS and mental/cog
threshold_list = ['0.5', '0.1', '0.05', '0.01', '0.001', '0.0001', '1e-05', '1e-06', '1e-07', '5e-08']
trait = 'cog'
N = 9
# trait = 'mental'
# N = 7
# 协变量
# gene_pc = pd.read_csv('D:/A/UK_Biobank/IDP/test2/all_IDP/results/GWAS/confound/gene_pc.csv')
gene_pc = pd.read_csv('D:/UKB/IDP/test2/all_IDP/results/GWAS/confound/gene_pc.csv')
# age_sex = pd.read_csv('D:/A/UK_Biobank/IDP/test2/fMRI/ICA/results/covariate/age_sex_502396.csv')
age_sex = pd.read_csv('age0_sex_502396.csv')
age_sex.columns = ['eid', 'age', 'sex']
cov = pd.merge(gene_pc, age_sex, on='eid')
cov_nonan = cov.dropna(axis=0, how='any')
cov_nonan['age2'] = cov_nonan['age'] ** 2
cov_nonan['age_sex'] = cov_nonan['age'] * cov_nonan['sex']
cov_nonan['age2_sex'] = cov_nonan['age'] ** 2 * cov_nonan['sex']
for threshold in threshold_list:
    print(threshold)
    my_prs1 = pd.read_csv('./correlation_' + trait + '/PRS_' + threshold + '/my_prs1.csv')
    my_prs2 = pd.read_csv('./correlation_' + trait + '/PRS_' + threshold + '/my_prs2.csv')
    data = pd.read_csv('./correlation_' + trait + '/' + trait + '_gene.csv')
    # ethnic = pd.read_csv('D:/A/UK_Biobank/IDP/test1/dMRI/genetic_ethnic/genetic_ethnic.csv', index_col=0)  # (502505, 2)
    ethnic = pd.read_csv('D:/UKB/IDP/test1/dMRI/genetic_ethnic/genetic_ethnic.csv', index_col=0) 
    my_ethnic = ethnic[ethnic['eid'].isin(my_prs1['old'])]
    white_British = my_ethnic.dropna(axis=0, how='any')
    my_prs1 = my_prs1[my_prs1['old'].isin(white_British['eid'])]
    my_prs2 = my_prs2[my_prs2['old'].isin(white_British['eid'])]
    data = data[data['eid'].isin(my_prs2['new'])]
    print(data.shape[0])

#     # 打印到文件‘result.txt’里
#     # with open('./correlation_' + trait + '/PRS_' + threshold + '/result1_white_outlier.txt', 'w') as f:
#     #     print('Pearson correlation:', file=f)
#     # with open('./correlation_' + trait + '/PRS_' + threshold + '/result_white_outlier_cov_regression.txt', 'w') as f:
#     with open('./correlation_' + trait + '/PRS_' + threshold + '/result_white_outlier_cov_regression_age0.txt', 'w') as f:  # age0
#     # with open('./correlation_' + trait + '/PRS_' + threshold + '/result_white_cov_regression.txt', 'w') as f:
#     # with open('./correlation_' + trait + '/PRS_' + threshold + '/result_white_cov1_regression.txt', 'w') as f:  # age, sex
#     # with open('./correlation_' + trait + '/PRS_' + threshold + '/result_white_cov2_regression.txt', 'w') as f:  # genepc20
#     # with open('./correlation_' + trait + '/PRS_' + threshold + '/result_white_cov3_regression.txt', 'w') as f:  # genepc10
#     # with open('./correlation_' + trait + '/PRS_' + threshold + '/result_white_outlier_cov3_regression.txt', 'w') as f:  # genepc10
#         for i in range(1, data.shape[1]):
#             # data删去离群值
#             data0 = copy.deepcopy(data)
#             delete_outlier(data0, i)
#             prs1 = my_prs1[my_prs1['new'].isin(data0['eid'])]
#             prs2 = my_prs2[my_prs2['new'].isin(data0['eid'])]
#             # prs删去离群值
#             delete_outlier(prs1, 2)
#             delete_outlier(prs2, 2)
#         #     data1 = data0[data0['eid'].isin(prs1['new'])]
#         #     data2 = data0[data0['eid'].isin(prs2['new'])]
#         #     print(data1.columns[i], file=f)
#         #     print(data1.shape[0], file=f)
#         #     print(data2.shape[0], file=f)
#         #     # 相关性
#         #     result1 = stats.pearsonr(prs1['score'], data1.iloc[:, i])
#         #     result2 = stats.pearsonr(prs2['score'], data2.iloc[:, i])
#         #     if result1[1] < 0.05/N:
#         #         print('mode1:', result1, file=f)
#         #     if result2[1] < 0.05/N:
#         #         print('mode2:', result2, file=f)
#         #     print('*' * 50, file=f)
#         # print('-' * 80, file=f)
#         #
#         # print('Spearman correlation:', file=f)
#         # for i in range(1, data.shape[1]):
#         #     result3 = stats.spearmanr(my_prs1['score'], data.iloc[:, i])
#         #     result4 = stats.spearmanr(my_prs2['score'], data.iloc[:, i])
#         #     print(data.columns[i], file=f)
#         #     if result3[1] < 0.05 / 9:
#         #         print('mode1:', result3, file=f)
#         #     if result4[1] < 0.05 / 9:
#         #         print('mode2:', result4, file=f)
#         #     print('*' * 50, file=f)
#
#             # 协变量
#             cov1 = cov_nonan[cov_nonan['eid'].isin(prs1['new'])]
#             cov2 = cov_nonan[cov_nonan['eid'].isin(prs2['new'])]
#             data1 = data0[data0['eid'].isin(cov1['eid'])]
#             data2 = data0[data0['eid'].isin(cov2['eid'])]
#             prs1 = prs1[prs1['new'].isin(cov1['eid'])]
#             prs2 = prs2[prs2['new'].isin(cov2['eid'])]
#             print(data1.columns[i], file=f)
#             print(data1.shape[0], file=f)
#             print(data2.shape[0], file=f)
#
#             # 回归分析
#             prs1.columns = ['eid', 'old', 'score']
#             prs2.columns = ['eid', 'old', 'score']
#             X1 = pd.merge(prs1, cov1, on='eid')
#             X2 = pd.merge(prs2, cov2, on='eid')
#             X1.drop(['eid', 'old'], axis=1, inplace=True)
#             X2.drop(['eid', 'old'], axis=1, inplace=True)
#             X1 = (X1 - X1.mean()) / X1.std()
#             X2 = (X2 - X2.mean()) / X2.std()
#             X1 = sm.add_constant(X1)
#             X2 = sm.add_constant(X2)
#             data1.set_index(X1.index, inplace=True)
#             data2.set_index(X2.index, inplace=True)
#             model1 = sm.OLS(data1.iloc[:, i], X1).fit()
#             model2 = sm.OLS(data2.iloc[:, i], X2).fit()
#             beta_values1 = model1.params
#             p_values1 = model1.pvalues
#             beta_values2 = model2.params
#             p_values2 = model2.pvalues
#             if p_values1['score'] < 0.05 / N:
#                 print('mode1_beta:', beta_values1['score'], file=f)
#                 print('mode1_p:', p_values1['score'], file=f)
#                 print('success!!!!')
#             if p_values2['score'] < 0.05 / N:
#                 print('mode2_beta:', beta_values2['score'], file=f)
#                 print('mode2_p:', p_values2['score'], file=f)
#                 print('success!!!!')
#             print('*' * 50, file=f)
#         print('-'*80, file=f)



threshold = '1e-05'
trait = 'cog'
i = 2
my_prs1 = pd.read_csv('./correlation_' + trait + '/PRS_' + threshold + '/my_prs1.csv')
my_prs2 = pd.read_csv('./correlation_' + trait + '/PRS_' + threshold + '/my_prs2.csv')
data = pd.read_csv('./correlation_' + trait + '/' + trait + '_gene.csv')
# 打印第i列的列名
print(data.columns[i])
ethnic = pd.read_csv('D:/UKB/IDP/test1/dMRI/genetic_ethnic/genetic_ethnic.csv', index_col=0)  # (502505, 2)
my_ethnic = ethnic[ethnic['eid'].isin(my_prs1['old'])]
white_British = my_ethnic.dropna(axis=0, how='any')
my_prs1 = my_prs1[my_prs1['old'].isin(white_British['eid'])]
my_prs2 = my_prs2[my_prs2['old'].isin(white_British['eid'])]
data = data[data['eid'].isin(my_prs2['new'])]
print(data.shape[0])
# data删去离群值
data0 = copy.deepcopy(data)
delete_outlier(data0, i)
print(data0.shape[0])

# prs1 = my_prs1[my_prs1['new'].isin(data0['eid'])]
prs2 = my_prs2[my_prs2['new'].isin(data0['eid'])]
# prs删去离群值
# delete_outlier(prs1, 2)
delete_outlier(prs2, 2)
# data1 = data0[data0['eid'].isin(prs1['new'])]
data2 = data0[data0['eid'].isin(prs2['new'])]
print(data2.shape[0])
# result1 = stats.pearsonr(prs1['score'], data1.iloc[:, i])
result2 = stats.pearsonr(prs2['score'], data2.iloc[:, i])
# result2 = stats.pearsonr(my_prs2['score'], data.iloc[:, i])
print(result2)
Mean = []
SD = []
N = []
# prs1['score_range_label'] = pd.cut(x=prs1['score'], bins=10, labels=range(10))
prs2['score_range_label'] = pd.cut(x=prs2['score'], bins=10, labels=range(10))
# prs2['score_range_label'].value_counts()
# my_prs2['score_range_label'] = pd.cut(x=my_prs2['score'], bins=10, labels=range(10))
for j in range(10):
    # my_prs = prs1[prs1['score_range_label'] == j]
    # my_data = data1[data1['eid'].isin(my_prs['new'])]
    my_prs = prs2[prs2['score_range_label'] == j]
    my_data = data2[data2['eid'].isin(my_prs['new'])]
    Mean.append(my_data[data.columns[i]].mean())
    SD.append(my_data[data.columns[i]].std())
    N.append(my_data.shape[0])
print(Mean)
print(SD)
print(N)
x = range(10)
plt.bar(x, Mean)
plt.show()

# 相关性分布图
sns.jointplot(data=data,
              x=my_prs1['score'],
              y=data['Numeric memory'],
              # y=data['Fluid intelligence'],
              kind="reg",
              height=5,
              color=(155/255, 62/255, 80/255))
# plt.savefig('./correlation_mental/PRS_' + threshold + '/mode1_Depression.eps', dpi=300)
plt.show()


# threshold = '0.5'
# trait = 'mental'
# i = 2
# my_prs1 = pd.read_csv('./correlation_' + trait + '/PRS_' + threshold + '/my_prs1.csv')
# my_prs2 = pd.read_csv('./correlation_' + trait + '/PRS_' + threshold + '/my_prs2.csv')
# data = pd.read_csv('./correlation_' + trait + '/' + trait + '_gene.csv')
# # 打印第i列的列名
# print(data.columns[i])
# ethnic = pd.read_csv('D:/UKB/IDP/test1/dMRI/genetic_ethnic/genetic_ethnic.csv', index_col=0)  # (502505, 2)
# my_ethnic = ethnic[ethnic['eid'].isin(my_prs1['old'])]
# white_British = my_ethnic.dropna(axis=0, how='any')
# my_prs1 = my_prs1[my_prs1['old'].isin(white_British['eid'])]
# my_prs2 = my_prs2[my_prs2['old'].isin(white_British['eid'])]
# data = data[data['eid'].isin(my_prs2['new'])]
# # data删去离群值
# data0 = copy.deepcopy(data)
# delete_outlier(data0, i)
#
# prs1 = my_prs1[my_prs1['new'].isin(data0['eid'])]
# # prs删去离群值
# delete_outlier(prs1, 2)
# data1 = data0[data0['eid'].isin(prs1['new'])]
# result1 = stats.pearsonr(prs1['score'], data1.iloc[:, i])
# # result2 = stats.pearsonr(my_prs2['score'], data.iloc[:, i])
# print(result1)
# Mean = []
# SD = []
# N = []
# prs1['score_range_label'] = pd.cut(x=prs1['score'], bins=10, labels=range(10))
# # prs2['score_range_label'].value_counts()
# # my_prs2['score_range_label'] = pd.cut(x=my_prs2['score'], bins=10, labels=range(10))
# for j in range(10):
#     my_prs = prs1[prs1['score_range_label'] == j]
#     my_data = data1[data1['eid'].isin(my_prs['new'])]
#     Mean.append(my_data[data.columns[i]].mean())
#     SD.append(my_data[data.columns[i]].std())
#     N.append(my_data.shape[0])
# print(Mean)
# print(SD)
# print(N)
# x = range(10)
# plt.bar(x, Mean)
# plt.show()
#
# # 相关性分布图
# sns.jointplot(data=data,
#               x=my_prs1['score'],
#               y=data['Trauma'],
#               # y=data['Depression'],
#               kind="reg",
#               height=5,
#               color=(155/255, 62/255, 80/255))
# # plt.savefig('./correlation_mental/PRS_' + threshold + '/mode1_Depression.eps', dpi=300)
# plt.show()


# CCA mode与认知的相关性
# IDP, gene, mental/cog的交集，中介分析

# # 中介分析
# data = pd.read_csv('./correlation_cog/cog_gene.csv')
# prs1 = pd.read_csv('./correlation_cog/my_prs1.csv')
# prs2 = pd.read_csv('./correlation_cog/my_prs2.csv')
# cca_mode = pd.read_csv('./mediation_analysis/CCA_U.csv')
# my_prs1 = prs1[prs1['new'].isin(cca_mode['eid'])]
# my_prs2 = prs2[prs2['new'].isin(cca_mode['eid'])]
# my_cca_mode = cca_mode[cca_mode['eid'].isin(my_prs2['new'])]
# my_cca_mode = my_cca_mode[['eid', 'mode1', 'mode2']]
# my_data = data[data['eid'].isin(my_cca_mode['eid'])]
# my_prs1.drop(['old'], axis=1, inplace=True)
# my_prs2.drop(['old'], axis=1, inplace=True)
# my_prs1.columns = ['eid', 'PRS1']
# my_prs2.columns = ['eid', 'PRS2']
# final_data = pd.merge(my_prs1, my_prs2, on='eid')
# final_data = pd.merge(final_data, my_cca_mode, on='eid')
# final_data = pd.merge(final_data, my_data, on='eid')
# print(final_data.shape[0])
# # confound = pd.read_csv('./mediation_analysis/Confounds_Subjects_19157.csv')
# # my_confound = confound[confound['eid'].isin(final_data['eid'])]
# # final_data = pd.merge(final_data, my_confound, on='eid')
# # final_data.to_csv('./mediation_analysis/mental/final_data.csv', index=False)
# confound = pd.read_csv(r'D:\A\UK_Biobank\IDP\test4\results\covariate\confound_19175.csv')
# my_confound = confound[confound['eid'].isin(final_data['eid'])]
# final_data = pd.merge(final_data, my_confound, on='eid')
# final_data.to_csv('./mediation_analysis/cog/final_data2.csv', index=False)

