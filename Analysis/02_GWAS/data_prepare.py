import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.model_selection import train_test_split


# Load cognition
# IDP_id = []
# with open("./behaviour/cognition/cognition_online.txt") as f:
#     for line in f:
#         if line != '\n':
#             line = line.strip()
#             IDP_id.append(line)  # 362

IDP_id = ['eid', '21003-0.0', '21003-2.0', '31-0.0']
# reader = pd.read_csv('E:/UK_Biobank/Download_669827/ukb669827.csv', iterator=True, chunksize=1000, encoding='gbk', low_memory=False)  # 按chunk读取UKB的IDP表,每1000条数据为一个chunk
reader = pd.read_csv('E:/UK_Biobank/original/ukb42416.csv', iterator=True, chunksize=1000, encoding='gbk', low_memory=False) 
df = None
n = 0
for chunk in reader:
    cnt_chunk = chunk[IDP_id]
    df = pd.concat([df, cnt_chunk], axis=0) if df is not None else cnt_chunk
    n += 1
    if n % 10 == 0:
        print(n)
# df.to_csv('D:/UKB/test4/results/covariate/new_confound.csv', index=False)
df.to_csv('D:/UKB/test4/results/covariate/old_confound.csv', index=False)


# # load IDP data
# data = pd.read_csv('D:/A/UK_Biobank/IDP/test2/all_IDP/results/IDP/norest_502396.csv')  # (502396, 363)
# IDP_id = []
# with open("./results/IDP/list.txt") as f:
#     for line in f:
#         if line != '\n':
#             line = line.strip()
#             IDP_id.append(line)
# idp = data[IDP_id]
# idp_nonan = idp.dropna(axis=0, how='any')
# idp_nonan.to_csv('./results/IDP/idp_sMRI.csv', index=False)  # 42301

# # mental health和IDP取交集
# idp = pd.read_csv('./results/IDP/idp_sMRI.csv')
# behaviour = pd.read_csv('D:/A/UK_Biobank/IDP/test3/mental_combine_2/results/feature/mental_zscore.csv')
# my_idp = idp[idp['eid'].isin(np.array(behaviour['eid']))]
# my_idp['eid'].to_csv('./results/eid/eid_overlap.csv', index=False)  # 19175

# # Prepare confounds
# eid = pd.read_csv('./results/eid/eid_overlap.csv')
# confound = pd.read_csv('D:/A/UK_Biobank/IDP/test2/fMRI/ICA/results/covariate/confound_502396.csv')
# my_confound = confound[confound['eid'].isin(np.array(eid['eid']))]
# my_confound1 = my_confound.drop(['25741-2.0', '25742-2.0'], axis=1)
# num = my_confound1.isna().sum()
# print(num)
# my_confound1.to_csv('./results/covariate/confound_' + str(my_confound1.shape[0]) + '.csv', index=False)  #19175

# # split exploratory sample(4/5) and confirmatory sample(1/5)
# data = pd.read_csv('./results/eid/eid_19175/eid_19175.csv')
# data_train, data_test = train_test_split(data, train_size=4/5, random_state=4)
# np.savetxt('./results/eid/eid_19175/random4/eid_train.csv', np.array(sorted(data_train['eid'])), delimiter=',')
# np.savetxt('./results/eid/eid_19175/random4/eid_test.csv', np.array(sorted(data_test['eid'])), delimiter=',')

# # eid mapping
# idp = pd.read_csv(r'D:\A\UK_Biobank\IDP\test4\results\IDP\idp_sMRI.csv')
# eid_mapping = pd.read_csv('D:/A/UK_Biobank/ID_mapping/eid_mapping.csv')
# my_eid = np.array(idp['eid'])
# eid = eid_mapping[eid_mapping['new'].isin(my_eid)]
# eid1 = eid.sort_values('new', inplace=False)
# eid1.to_csv('./results/GWAS/eid/eid_mapping.csv', index=False)

# # 去掉ethnic的nan，并找出white british被试
# ethnic = pd.read_csv('D:/A/UK_Biobank/IDP/test1/dMRI/genetic_ethnic/genetic_ethnic.csv', index_col=0)  # (502505, 2)
# eid = pd.read_csv('./results/GWAS/eid/eid_nosex.csv', header=None)
# my_ethnic = ethnic[ethnic['eid'].isin(np.array(eid[0]))]
# white_British = my_ethnic.dropna(axis=0, how='any')
# white_British['eid'].to_csv('./results/GWAS/eid/eid_white.csv', index=False)

# # 找unrelated subjects: 0 No kinship found
# kinship = pd.read_csv('D:/A/UK_Biobank/IDP/test1/dMRI/genetic_kinship/genetic_kinship.csv', index_col=0)  # (502505, 2)
# eid = pd.read_csv('./results/GWAS/eid/eid_white.csv')
# my_kinship = kinship[kinship['eid'].isin(np.array(eid['eid']))]
# no_kinship = my_kinship[my_kinship['22021-0.0'] == 0]
# no_kinship['eid'].to_csv('./results/GWAS/eid/eid_nokinshp.csv', index=False)

# # largest unrelated subset
# ID = []
# with open("D:/A/UK_Biobank/IDP/test2/all_IDP/results/GWAS/kinship/ukb62352_rel_chr10_s488264.dat") as f:
#     for line in f:
#         line = line.strip().split()
#         ID.append(line[0:2])
# ID = ID[1:]
# eid = pd.read_csv('./results/GWAS/eid/eid_white.csv')
# eid = np.array(eid)
# # 如果连接一条边的两个顶点中都属于我的样本，则称这条边为 E2
# my_ID = [i for i, a in enumerate(ID) if int(a[0]) in eid and int(a[1]) in eid]
# ID = np.array(ID)
# my_subjects = ID[my_ID]
# np.savetxt('./results/GWAS/kinship/my_pairwise_and/my_pairwise_' + str(my_subjects.shape[0]) + '.csv',
#            my_subjects, delimiter=',', fmt='%s')
# # 如果连接一条边的两个顶点中只有一个点属于我的样本，则称这条边为 E1
# my_ID2 = [i for i, a in enumerate(ID) if sum([int(a[0]) in eid, int(a[1]) in eid]) == 1]
# my_subjects2 = ID[my_ID2]
# np.savetxt('./results/GWAS/kinship/my_pairwise_or-and/my_pairwise_or-and_' + str(my_subjects2.shape[0])
#            + '.csv', my_subjects2, delimiter=',', fmt='%s')
# # E1去掉重复的点
# ID_or_and = my_subjects2
# my_ID_1 = [a[0] for a in ID_or_and if int(a[0]) in eid]
# my_ID_2 = [a[1] for a in ID_or_and if int(a[1]) in eid]
# my_ID_all = set(my_ID_1 + my_ID_2)
# np.savetxt('./results/GWAS/kinship/my_pairwise_or-and/no_kinship_' + str(len(my_ID_all)) + '.csv',
#            np.array(list(my_ID_all)), delimiter=',', fmt='%s')

# # 按照边列表的形式读入文件，生成无向图
# list = pd.read_csv('./results/GWAS/kinship/my_pairwise_and/my_pairwise_1001.csv', header=None)
# list1 = np.array(list)
# list2 = [[str(x) for x in row] for row in list1]
# n = 0
# name2id = {}
# id2name = {}
# list3 = []
# for row in list2:
#     tmp = []
#     tmp1 = []
#     for x in row:
#         if x not in name2id.keys():
#             name2id[x] = n
#             id2name[n] = x
#             n += 1
#         tmp.append(name2id[x])
#     list3.append(tmp)
# np.savetxt('./results/GWAS/kinship/my_pairwise_and/my_pairwise_1001_1.txt', np.array(list3), fmt='%d')

# import igraph as ig

# g = ig.Graph.Read_Edgelist("./results/GWAS/kinship/my_pairwise_and/my_pairwise_1001_1.txt",
#                            directed=False)
# layout = g.layout("tree")
# image = ig.plot(g, layout=layout)
# image.save('./results/GWAS/kinship/my_pairwise_and/igraph.png')

# # 用 networkx 找 maximal independent set
# import networkx as nx

# g = nx.Graph()
# ID = []
# with open("./results/GWAS/kinship/my_pairwise_and/my_pairwise_1001_1.txt") as f:
#     for line in f:
#         line = line.strip().split()
#         line1 = [int(x) for x in line]
#         ID.append(tuple(line1))
# g.add_edges_from(ID)
# eid = pd.read_csv('./results/GWAS/eid/eid_white.csv')
# eid = np.array(eid)
# all_id = []
# for j in range(10000):
#     print(j)
#     max_subset = nx.maximal_independent_set(g)
#     print(len(max_subset))
#     max_id = [id2name[max_subset[i]] for i in range(len(max_subset))]
#     my_id = [max_id[i] for i in range(len(max_id)) if int(max_id[i]) in eid]
#     all_id.append(my_id)
# max_len = len(max(all_id, key=len))
# print(max_len)
# final_id = max(all_id, key=len)
# final_id1 = [int(x) for x in final_id]
# np.savetxt('./results/GWAS/kinship/my_pairwise_and/no_kinship_' + str(len(final_id)) + '.csv',
#            np.array(final_id1), delimiter=',')

# # E1和E2中有重复的样本
# ID1 = pd.read_csv('./results/GWAS/kinship/my_pairwise_and/no_kinship_951.csv', header=None)
# ID2 = pd.read_csv('./results/GWAS/kinship/my_pairwise_or-and/no_kinship_10350.csv', header=None)
# ID3 = pd.read_csv('./results/GWAS/eid/eid_nokinshp.csv')
# x1 = [int(a) for a in np.array(ID1)]
# x2 = [int(a) for a in np.array(ID2)]
# x3 = [int(a) for a in np.array(ID3['eid'])]
# y = x1 + x2
# y_clean = set(x1 + x2)
# # 去掉ID1和ID2里有亲缘关系的ID
# pairwise1 = pd.read_csv('./results/GWAS/kinship/my_pairwise_and/my_pairwise_1001.csv', header=None)
# pairwise1 = np.array(pairwise1)
# pair_ID_1 = [a[1] for i, a in enumerate(pairwise1) if a[0] in x1]  # 323
# pair_ID_2 = [a[0] for i, a in enumerate(pairwise1) if a[1] in x1]  # 334
# pair_ID = pair_ID_1 + pair_ID_2  # 657, 与647顶点相连的点有 657=323+334个,里面有重复的点
# x4 = [a for a in x2 if not (a in pair_ID)]  # 8235, 8516里去掉与647顶点相连的点（657），还剩8235个点
# np.savetxt('./results/GWAS/kinship/my_pairwise_or-and/no_kinship_' + str(len(x4)) + '.csv',
#            np.array(x4), delimiter=',', fmt='%s')
# y1 = x1 + x3 + x4  # 28512
# y1_final = set(x1 + x3 + x4)  # 28235
# np.savetxt('./results/GWAS/kinship/eid_nokinship_' + str(len(y1_final)) + '.csv',
#            np.array(list(y1_final)), delimiter=',', fmt='%s')


# # GWAS: prepare confound
# # eid = pd.read_csv('D:/A/UK_Biobank/ID_mapping/eid_mapping.csv')
# # eid_nokinship = pd.read_csv('./results/GWAS/kinship/eid_nokinship_34869.csv', header=None)
# # eid_gwas = eid[eid['old'].isin(np.array(eid_nokinship[0]))]
# # eid_gwas = eid_gwas.sort_values(by='new')
# # eid_gwas.to_csv('./results/GWAS/eid/eid_gwas_new_old.csv', index=False)
# # eid_gwas['new'].to_csv('./results/GWAS/eid/eid_gwas_new.csv', index=False, header=None)
# eid = pd.read_csv('./results/GWAS/eid/eid_gwas_new.csv', header=None)
# confound = pd.read_csv('D:/A/UK_Biobank/IDP/test2/fMRI/ICA/results/covariate/confound_502396.csv')
# confound1 = confound[confound['eid'].isin(np.array(eid[0]))]
# my_confound1 = confound1.drop(['20400-0.0', '25741-2.0', '25742-2.0'], axis=1)
# num = my_confound1.isna().sum()
# print(num)
# my_confound1.to_csv('./results/GWAS/covariate/confound_' + str(confound1.shape[0]) + '.csv', index=False)
# gene_pc = pd.read_csv('D:/A/UK_Biobank/IDP/test2/all_IDP/results/GWAS/confound/gene_pc.csv')  # (502396, 41)
# gene_pc1 = gene_pc[gene_pc['eid'].isin(np.array(eid[0]))]
# gene_pc1.to_csv('./results/GWAS/covariate/gene_pc.csv', index=False)

# # prepare GWAS data
# pheno = pd.read_csv(r'.\results\GWAS\mental\GWAS_behaviour1.csv', header=None)
# eid = pd.read_csv(r'D:\A\UK_Biobank\IDP\test3\mental_combine_2\behaviour\results\GWAS\eid\eid_gwas_new_old.csv')
# pheno.columns = ['feature1']
# pheno.insert(0, 'EID', eid['old'].tolist())
# pheno.insert(1, 'IID', eid['old'].tolist())
# s = pheno['feature1']
# # 离群值
# threshold = s.std() * 3
# outliers = []
# for y in s:
#     if np.abs(y - s.mean()) > threshold:
#         outliers.append(y)
#         index = pheno[pheno.feature1 == y].index.to_list()
#         pheno.drop(index, inplace=True)
# print(len(outliers))
# # pheno.to_csv('./results/GWAS/mental/pheno2.txt', sep='\t', header=True, index=False)
# # pheno[['EID', 'IID']].to_csv('./results/GWAS/mental/eid2_' + str(pheno.shape[0]) + '.txt',
# #                             sep='\t', header=True, index=False)

# # for i in range(2, 23):
# #     a = '/mnt/data2/home02/khu/UKB_GWAS/CCA_combine/mental_sMRI/GWAS_mental/snp_filter_missnp/ukb22828_c' + str(i) + '_4'
# #     print(a + '.bed' + ' ' + a + '.bim' + ' ' + a + '.fam')
