import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from lifelines import KaplanMeierFitter
from lifelines.utils import median_survival_times
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from lifelines import WeibullFitter, ExponentialFitter, LogNormalFitter, LogLogisticFitter

# prepare data
# data = []
# with open('field.txt') as f:
#     for line in f:
#         if line != '\n':
#             line = line.strip()
#             data.append(line)
field = ['eid', '53-0.0', '21003-0.0']
reader = pd.read_csv('D:/A/UK_Biobank/Download_669827/ukb669827.csv', iterator=True, chunksize=1000, encoding='gbk', low_memory=False)  # 按chunk读取UKB的IDP表,每1000条数据为一个chunk
df = None
n = 0
for chunk in reader:
    cnt_chunk = chunk[field]
    df = pd.concat([df, cnt_chunk], axis=0) if df is not None else cnt_chunk
    n += 1
    if n % 10 == 0:
        print(n)
df.to_csv('./date_age0.csv', index=False)

# prepare my data
data = pd.read_csv('data.csv')
date = pd.read_csv('date0.csv')
data = pd.merge(data, date, on='eid')
cca_mode = pd.read_csv('GWAS_IDP.csv')
my_data = data[data['eid'].isin(cca_mode['new'])]
cca_mode = cca_mode[cca_mode['new'].isin(my_data['eid'])]
my_data = my_data.drop(['40000-1.0'], axis=1)
print(my_data.isnull().sum())
data_nonan = my_data.dropna(axis=0, how='any')
my_data['40000-0.0'] = my_data['40000-0.0'].fillna('2021/12/1')
my_data['40023-0.0'] = my_data['40023-0.0'].fillna(0)
my_data['40000-0.0'] = pd.to_datetime(my_data['40000-0.0'], format='%Y/%m/%d')
# my_data['53-2.0'] = pd.to_datetime(my_data['53-2.0'], format='%Y/%m/%d')
# my_data['time'] = (my_data['40000-0.0'] - my_data['53-2.0']).dt.days
my_data['53-0.0'] = pd.to_datetime(my_data['53-0.0'], format='%Y/%m/%d')
my_data['time'] = (my_data['40000-0.0'] - my_data['53-0.0']).dt.days
my_data.columns = ['new', 'sex', 'age', 'date_begin2', 'date_died', 'status', 'date_begin0','time']
final_data = pd.merge(cca_mode, my_data, on='new')
final_data.to_csv('my_data_date0.csv', index=False)


# # prepare original mental health data
# my_data = pd.read_csv('my_data.csv')
# mental = pd.read_csv(r'D:\A\UK_Biobank\IDP\test3\mental_combine_2\results\feature\mental_zscore.csv')
# my_mental = mental[mental['eid'].isin(my_data['new'])]
# my_data2 = my_data[my_data['new'].isin(my_mental['eid'])]
# my_mental.columns = ['new', 'Anxiety', 'Trauma', 'Depression', 'Psychotic_experience', 'Self_harm', 'Mania',
#                      'Mental_distress']
# data = pd.merge(my_data2, my_mental, on='new')
# confound = pd.read_csv('D:/A/UK_Biobank/IDP/test2/fMRI/ICA/results/covariate/confound_502396.csv')
# mental_date = confound[['eid', '20400-0.0']]
# mental_date.columns = ['new', 'mental_date']
# data = pd.merge(data, mental_date, on='new')
# data['date_died'] = pd.to_datetime(data['date_died'], format='%Y/%m/%d')
# data['mental_date'] = pd.to_datetime(data['mental_date'], format='%Y/%m/%d')
# data['mental_time'] = (data['date_died'] - data['mental_date']).dt.days
# data.to_csv('my_data2.csv', index=False)

# # prepare PRS data
# # eid_rest_gene_white = pd.read_csv('../post_analysis/rest_gene_eid_white_British.csv')
# # data = pd.read_csv('data.csv')
# # date = pd.read_csv('date_age0.csv')
# # data = pd.merge(data, date, on='eid')
# # my_data = data[data['eid'].isin(eid_rest_gene_white['new'])]
# # my_data = my_data.drop(['40000-1.0'], axis=1)
# # print(my_data.isnull().sum())
# # my_data['40000-0.0'] = my_data['40000-0.0'].fillna('2021/12/1')
# # my_data['40023-0.0'] = my_data['40023-0.0'].fillna(0)
# # my_data['40000-0.0'] = pd.to_datetime(my_data['40000-0.0'], format='%Y/%m/%d')
# # # my_data['53-2.0'] = pd.to_datetime(my_data['53-2.0'], format='%Y/%m/%d')
# # # my_data['time'] = (my_data['40000-0.0'] - my_data['53-2.0']).dt.days
# # my_data['53-0.0'] = pd.to_datetime(my_data['53-0.0'], format='%Y/%m/%d')
# # my_data['time'] = (my_data['40000-0.0'] - my_data['53-0.0']).dt.days
# # my_data.columns = ['new', 'sex', 'age', 'date_begin2', 'date_died', 'status', 'date_begin0', 'age0', 'time']
# # my_data = pd.merge(eid_rest_gene_white, my_data, on='new')
# # my_data.to_csv('data_time0.csv', index=False)
# my_data = pd.read_csv('data_time0.csv')
# threshold = '5e-08'
# prs1 = pd.read_csv('../post_analysis/PRS_plink/mode1/PRS_' + threshold + '.csv')
# prs2 = pd.read_csv('../post_analysis/PRS_plink/mode2/PRS_' + threshold + '.csv')
# prs1.columns = ['old', 'score1']
# prs2.columns = ['old', 'score2']
# final_data = pd.merge(my_data, prs1, on='old')
# final_data = pd.merge(final_data, prs2, on='old')
# final_data.to_csv(threshold + '_PRS_date0.csv', index=False)




# survival analysis
# threshold = '5e-08'
# data = pd.read_csv(threshold + '_PRS_date0.csv')
# data = data[['status', 'time', 'score1', 'score2']]
# # zscore标准化
# data['score1'] = (data['score1'] - data['scmode_date2ore1'].mean()) / data['score1'].std()
# data['score2'] = (data['score2'] - data['score2'].mean()) / data['score2'].std()

data = pd.read_csv('mode_date2.csv')
data['time_years'] = data['time'] / 365
data.to_csv('mode_date2_years.csv')
data = data[['mode1', 'mode2', 'sex', 'age', 'status', 'time_years']]  # age, sex, mode
# data = data[['mode1', 'mode2', 'status', 'time']]  # mode
# data = data[['sex', 'status', 'time']]  # mode
# data = pd.read_csv('my_data2.csv')
# data = data[['mode1', 'mode2', 'sex', 'age', 'status', 'mental_time', 'Anxiety', 'Trauma', 'Depression',
#              'Psychotic', 'Self_harm', 'Mania', 'Mental_distress']]  # mental health, mode
# data = data[['sex', 'age', 'status', 'mental_time', 'Anxiety', 'Trauma', 'Depression',
#              'Psychotic', 'Self_harm', 'Mania', 'Mental_distress']]  # mental health
# print(min(data['score1']))
# print(max(data['score1']))
# print(min(data['score2']))
# print(max(data['score2']))




# Cox Proportional Hazard Model (Semi-Parametric)
cph = CoxPHFitter()
cph.fit(data, duration_col='time_years', event_col='status')
# # cph.fit(data, duration_col='mental_time', event_col='status')
# cph.print_summary()
# plt.subplots(figsize=(10, 6))
# cph.plot()
# plt.savefig('./results/age_sex_mode/cph', format='eps', dpi=600)
# plt.show()

# # Plot Partial Effects on Outcome (Cox-PH Regression)
cph.plot_partial_effects_on_outcome(['score1', 'score2'], values=[
    [3, 0],
    [1, 0],
    [-1, 0],
    [-3, 0],
    [0, -3],
    [0, -1],
    [0, 1],
    [0, 3]
], cmap='coolwarm')
# plt.savefig('./results/age0_sex_PRS/effect2_' + threshold + '.png')
# plt.show()

# cph.plot_partial_effects_on_outcome(covariates='mode1', values=[-5, -3, -1, 1, 3, 5], cmap='Blues', )
# plt.ylim((0.95, 1))
# plt.savefig('./results/mode1.png')
# plt.show()
# cph.plot_partial_effects_on_outcome(covariates='mode2', values=[-5, -3, -1, 1, 3, 5], cmap='Reds')
# plt.ylim((0.95, 1))
# plt.savefig('./results/mode2.png')
# plt.show()

# Plot Partial Effects on Outcome (Cox-PH Regression)
cph.plot_partial_effects_on_outcome(['mode1', 'mode2'], values=[
    [5, 0],
    # [3, 0],
    # [1, 0],
    # [-1, 0],
    # [-3, 0],
    # [-5, 0],
    # [0, -5],
    # [0, -3],
    # [0, -1],
    # [0, 1],
    # [0, 3],
    [0, 5]
], cmap='coolwarm_r')

# 设置标签
plt.xlabel('Time (years)')
plt.ylabel('Survival probability')
# 调节文字大小
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
# 调节图例大小
plt.legend(fontsize=8)
plt.savefig('./results/age_sex_mode/mode_years_eg.svg', dpi=600)
plt.show()



# # Check Proportional Hazard Assumption
# cph.check_assumptions(data, p_value_threshold=0.05)
results = proportional_hazard_test(cph, data, time_transform='rank')
# results.print_summary(decimals=3, model="untransformed variables")

# # correlation: PRS and mode
# mode = pd.read_csv(r'D:\A\UK_Biobank\IDP\test4\post_analysis\PRSice\white\res_cov_top40pc\GWAS_IDP\GWAS_IDP_white.csv')
# threshold = '5e-08'
# prs1 = pd.read_csv('../post_analysis/PRS_plink/mode1/PRS_' + threshold + '.csv')
# prs2 = pd.read_csv('../post_analysis/PRS_plink/mode2/PRS_' + threshold + '.csv')
# prs1.columns = ['old', 'score1']
# prs2.columns = ['old', 'score2']
# prs1 = prs1[prs1['old'].isin(mode['old'])]
# prs2 = prs2[prs2['old'].isin(mode['old'])]
# final_data = pd.merge(mode, prs1, on='old')
# final_data = pd.merge(final_data, prs2, on='old')
# result1 = stats.spearmanr(final_data['mode1'], final_data['score1'])
# result2 = stats.spearmanr(final_data['mode2'], final_data['score1'])
