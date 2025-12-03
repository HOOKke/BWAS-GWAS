clear all; close all; clc

addpath('/mnt/data0/home/khu/CCA/MatlabScripts/FSLNets')
addpath('/mnt/data0/home/khu/CCA/MatlabScripts/FSL')
addpath('/mnt/data0/home/khu/CCA/MatlabScripts/')
% INPUT = 'eid_19157.csv';
INPUT = 'eid_test.csv';
% INPUT = 'eid_overlap.csv';


%% Load confound
load(sprintf('%s/Confounds_Subjects_%s.mat','/mnt/data0/home/khu/CCA/Confound/confound_test3/mental_sMRI/eid_19157/random4',INPUT(5:end-4)));
for i = 1:size(conf,2)
    if any(conf(:,i)) == 1
        conf(:,i) = nets_normalise(conf(:,i));
    end
end
Pconf = pinv(conf);

%% IDPs and behaviour
fprintf('Processing IDPs\n');
load('./eid_19157/random4/PC_70/new_results2/Eigs_train');
S = load(sprintf('%s/%s','/mnt/data0/home/khu/CCA/Data/mental_sMRI/eid_19157/random4',INPUT));
% S = importdata(sprintf('%s/%s','/mnt/data0/home/khu/CCA/Data/mental_sMRI/eid_19157',INPUT));
% S = S.data;
IDP = importdata('/mnt/data0/home/khu/CCA/Data/mental_sMRI/idp_sMRI.csv');
IDP = IDP.data;
subs = IDP(:,1); IDP = IDP(:,2:end);
[~,s,~] = intersect(subs,S); IDP = IDP(s,:);
IDP = nets_normalise(IDP);
IDP = nets_normalise(IDP-conf*(Pconf*IDP));
SCORE_PROJECT = IDP*REST_COEFF;
REST_OUT = SCORE_PROJECT(:,1:Nkeep_rest);  clear SCORE_PROJECT

behaviour = importdata('/mnt/data0/home/khu/CCA/Data/mental_sMRI/mental_zscore.csv');
behaviour = behaviour.data;
subs = behaviour(:,1); behaviour = behaviour(:,2:end);
[~,s,~] = intersect(subs,S); behaviour = behaviour(s,:);
behaviour = nets_normalise(behaviour);
behaviour = nets_normalise(behaviour-conf*(Pconf*behaviour));

%% Save results
sprintf('Saving results\n');
save(sprintf('./eid_19157/random4/PC_70/new_results2/CCA_inputs_brain_%s.mat',INPUT(5:end-4)),'REST_OUT', 'IDP');
save(sprintf('./eid_19157/random4/PC_70/new_results2/CCA_inputs_behaviour_%s.mat',INPUT(5:end-4)),'behaviour');
dlmwrite(sprintf('./eid_19157/random4/PC_70/new_results2/CCA_inputs_brain_%s.csv',INPUT(5:end-4)),IDP,'precision',25);

