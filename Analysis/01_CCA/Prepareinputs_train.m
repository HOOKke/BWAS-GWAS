clear all; close all; clc

addpath('/mnt/data0/home/khu/CCA/MatlabScripts/FSLNets')
addpath('/mnt/data0/home/khu/CCA/MatlabScripts/FSL')
addpath('/mnt/data0/home/khu/CCA/MatlabScripts/')
INPUT = 'eid_train.csv';
Checks = zeros(1,3); 
Nkeep = 70; % keep number of PCs that explain at least 50% of variance

%% Load confound
load(sprintf('%s/Confounds_Subjects_%s.mat','/mnt/data0/home/khu/CCA/Confound/confound_test3/mental_sMRI/eid_19157/random4/new_15323',INPUT(5:end-4)));
for i = 1:size(conf,2)
    if any(conf(:,i)) == 1
        conf(:,i) = nets_normalise(conf(:,i));
    end
end
Pconf = pinv(conf);

%% IDPs and behaviour
fprintf('Processing IDPs\n');
S = load(sprintf('%s/%s','/mnt/data0/home/khu/CCA/Data/mental_sMRI/eid_19157/random4',INPUT));
S(8608) = [];
S(10444) = [];
NkeepMax = floor(length(S)/10/2);

IDP = importdata('/mnt/data0/home/khu/CCA/Data/mental_sMRI/idp_sMRI.csv');
IDP = IDP.data;
subs = IDP(:,1); IDP = IDP(:,2:end);
[~,s,~] = intersect(subs,S); IDP = IDP(s,:);
IDP = nets_normalise(IDP);
IDP = nets_normalise(IDP-conf*(Pconf*IDP));
[REST_COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(IDP);
% % 绘制各成分累计解释的方差图
% plot(cumsum(EXPLAINED),'ro-');
% xlabel('主成分数量');
% ylabel('累计方差解释比例');
Vexp = cumsum(EXPLAINED); 
% Nkeep = find(Vexp>=variance_threshold,1,'first');
% if Nkeep > NkeepMax; Nkeep = NkeepMax; end
REST_OUT = SCORE(:,1:Nkeep);
Checks(1,1) = size(IDP,2); Checks(1,2) = Nkeep; Checks(1,3) = Vexp(Nkeep);
Nkeep_rest = Nkeep;
clear COEFF LATENT TSQUARED EXPLAINED MU Nkeep

behaviour = importdata('/mnt/data0/home/khu/CCA/Data/mental_sMRI/mental_zscore.csv');
behaviour = behaviour.data;
subs = behaviour(:,1); behaviour = behaviour(:,2:end);
[~,s,~] = intersect(subs,S); behaviour = behaviour(s,:);
behaviour = nets_normalise(behaviour);
behaviour = nets_normalise(behaviour-conf*(Pconf*behaviour));

behaviour = importdata('/mnt/data0/home/khu/CCA/Data/mental_sMRI/cognition_2_nonan_norm.csv');
behaviour = behaviour.data;
subs = behaviour(:,1); behaviour = behaviour(:,2:end);
[~,s,~] = intersect(subs,S); behaviour = behaviour(s,:);
CCA = load('./eid_19157/random4/PC_70/new_results2/CCA_results_train.mat');
U = CCA.U;
[~,s2,~] = intersect(S, subs); idp_mode = U(s2,:);

%% Save results
sprintf('Saving results\n');
save(sprintf('/mnt/data0/home/khu/CCA/CCA/CCA_test3/eid_19157/random4/PC_70/new_results2/CCA_inputs_brain_%s.mat',INPUT(5:end-4)),'REST_OUT', 'IDP', 'Checks');
save(sprintf('/mnt/data0/home/khu/CCA/CCA/CCA_test3/eid_19157/random4/PC_70/new_results2/Eigs_%s.mat',INPUT(5:end-4)),'REST_COEFF','Nkeep_rest');
save(sprintf('/mnt/data0/home/khu/CCA/CCA/CCA_test3/eid_19157/random4/PC_70/new_results2/CCA_inputs_behaviour_%s.mat',INPUT(5:end-4)),'behaviour');
