clear all; close all; clc

addpath('/mnt/data0/home/khu/CCA/MatlabScripts/')

% Run CCA on all subjects
fprintf('Running CCA on all subjects \n');
load('./eid_19157/random4/PC_70/new_results2/CCA_inputs_behaviour_train.mat');
load('./eid_19157/random4/PC_70/new_results2/CCA_inputs_brain_train.mat');

[pfwer,r,A,B,U,V] = permcca(REST_OUT,behaviour,1000);
save('./eid_19157/random4/PC_70/new_results2/CCA_results_train.mat','pfwer','r','A','B','U','V'); 