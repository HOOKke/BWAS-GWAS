clear all; close all; clc

addpath('/mnt/data0/home/khu/CCA/MatlabScripts/FSLNets')
addpath('/mnt/data0/home/khu/CCA/MatlabScripts/FSL')
addpath('/mnt/data0/home/khu/CCA/MatlabScripts/')
CV = 1;


%% Load data
CCA = load('./eid_19157/random4/PC_70/new_results2/CCA_results_train.mat');
load('./eid_19157/random4/PC_70/new_results2/CCA_inputs_brain_test.mat');
load('./eid_19157/random4/PC_70/new_results2/CCA_inputs_behaviour_test.mat');
% load('./eid_19157/random4/PC_70/CCA_inputs_brain_19157.mat');
% load('./eid_19157/random4/PC_70/CCA_inputs_behaviour_19157.mat');
Urepl = REST_OUT * CCA.A;
Vrepl = behaviour * CCA.B;
% save('./eid_19157/random4/PC_70/CCA_U.mat', 'Urepl');
% save('./eid_19157/random4/PC_70/CCA_V.mat', 'Vrepl');
% dlmwrite('./eid_19157/random4/PC_70/CCA_U.csv',Urepl,'precision',25);
% dlmwrite('./eid_19157/random4/PC_70/CCA_V.csv',Vrepl,'precision',25);
for i = 1:size(Urepl,2)
    [r(i),p(i)] = corr(Urepl(:,i),Vrepl(:,i));
end
p_FDR = mafdr(p,'BHFDR', true);
p_fwer = p*length(p);
save('./eid_19157/random4/PC_70/new_results2/CCA_results_test.mat','r','p','p_FDR','p_fwer');             


%% correlation with cognition
load('/mnt/data0/home/khu/CCA/CCA/CCA_test3/eid_19157/random4/PC_70/Cognition.mat')
load('/mnt/data0/home/khu/CCA/CCA/CCA_test3/eid_19157/random4/PC_70/Cognition_IDP_mode.mat')
[r1,p1] = corr(behaviour,idp_mode(:,1));
[r2,p2] = corr(behaviour,idp_mode(:,2));

%% Correlation with eigenvectors rather than IDPs (for reviewer)
N = size(IDP, 2);
load('./eid_19157/random4/PC_70/new_results2/CCA_inputs_brain_train.mat');
[r1,p1] = corr(CCA.U(:,CV), IDP);
load('./eid_19157/random4/PC_70/new_results2/CCA_inputs_brain_test.mat');
[r2,p2] = corr(Urepl(:,CV), IDP);
ind = p1<0.05/N/7 & p2<0.05/N/7;
TEST(:,1) = r1(ind)';
TEST(:,2) = r2(ind)';
% TEST = abs(TEST);

load('IDP_name.mat');
for i = 1:length(ind)
    if ind(i)==0
        IDP_name{i} = {};
    end
end
IDP_name(cellfun(@isempty,IDP_name)) = [];

[TEST0, sort_ind] = sortrows(TEST,'ascend');
name = cell(1,length(IDP_name));
for i = 1:length(IDP_name)
    name{i} = IDP_name{sort_ind(i)};
end
name=name';
save(['./eid_19157/random4/PC_70/loadings/IDP_loading', num2str(CV), '.mat'],'TEST0','name');
% save(['./loadings/IDP_loading', num2str(CV), '(abs).mat'],'TEST0','name');
% icc_coef = ICC(zscore(TEST0),'1-1');
% save(['./loadings/ICC/IDP_ICC', num2str(CV), '.mat'],'icc_coef');


% figure; set(gcf,'Position',[100 100 1300 950],'PaperPositionMode','auto')
% barh(TEST0(:,1),0.5,'FaceColor','r')
% hold on
% barh(TEST0(:,2),0.25,'FaceColor','b')
% legend({'Correlation with UV in exploratory CCA sample','Correlation with UV in replication sample'},'Location','SouthEast')
% xlabel('Cross-subject correlation with CCA score (UV)'); 
% title('IDP CCA correlations');
% set(gca, 'ytick',1:1:length(IDP_name), 'YTickLabel', name, 'FontSize',4);
% print(gcf,['./test_split_5_70/loadings/IDP_loading', num2str(CV)],'-dpng','-r300');



%% Behavioral correlations
load('./eid_19157/random4/PC_60/CCA_inputs_behaviour_train.mat');
R = corr(behaviour,CCA.V(:,CV));
load('./eid_19157/random4/PC_60/CCA_inputs_behaviour_test.mat');
R2 = corr(behaviour,Vrepl(:,CV));
R_all = [R, R2];
save(['./eid_19157/random4/PC_60/loadings/behaviour_loading', num2str(CV), '.mat'],'R_all');
% Anxiety	Trauma	Depression	Psychotic_experience	Self_harm	Mania	Mental_distress

% R = abs(R);
% R2 = abs(R2);
% save(['./loadings/behaviour_loading', num2str(CV), '(abs).mat'],'R','R2');
% icc_coef = ICC(zscore(R_all),'1-1');
% save(['./loadings/ICC/behaviour_ICC', num2str(CV), '.mat'],'icc_coef');


% figure
% % spider_plot([abs(R2) abs(R)]',...
% spider_plot([R2 R]',...
%     'AxesLabels', {'Wellbeing','Anxiety','Trauma','Depression','Psychotic experience','Self harm','Mania','Mental distress'},...
%     'AxesInterval', 4,...
%     'AxesPrecision', 2,...
%     'AxesLimits', [-1,-1,-1,-1,-1,-1,-1,-1;1,1,1,1,1,1,1,1],...
%     'Color', [0, 0, 1; 1, 0, 0],...
%     'AxesFontSize', 8,...
%     'LabelFontSize', 10);
% legend({'Correlation with UV in replication sample','Correlation with UV in exploratory CCA sample'},'Location','SouthEast')
% title('behaviour CCA correlations');
% % print(gcf,['./loadings/behaviour_loading', num2str(CV), '(abs)'],'-dpng','-r300');
% print(gcf,['./loadings/behaviour_loading', num2str(CV)],'-dpng','-r300');


% name0 = {'Anxiety','Trauma','Depression','Psychotic experience','Self harm','Mania','Mental distress'};
% [R_all0, sort_ind] = sortrows(R_all,'ascend');
% name1 = cell(1,length(name0));
% for i = 1:length(name0)
%     name1{i} = name0{sort_ind(i)};
% end
% % a = name0';
% figure; set(gcf,'Position',[100 100 1300 950],'PaperPositionMode','auto')
% barh(R_all0(:,1),0.5,'FaceColor','r')
% hold on
% barh(R_all0(:,2),0.25,'FaceColor','b')
% legend({'Correlation with UV in exploratory CCA sample','Correlation with UV in replication sample'},'Location','SouthEast')
% xlabel('Cross-subject correlation with CCA score (UV)'); 
% title('behaviour CCA correlations');
% set(gca, 'ytick',1:1:length(name1), 'YTickLabel', name1, 'FontSize',15);
% print(gcf,['./test_split_5_70/loadings/behaviour_loading', num2str(CV)],'-dpng','-r300');


% R_all = -R_all;
% name0 = {'Anxiety','Trauma','Depression','Psychotic experience','Self harm','Mania','Mental distress'};
% [R_all0, sort_ind] = sortrows(R_all,'descend');
% name1 = cell(1,length(name0));
% for i = 1:length(name0)
%     name1{i} = name0{sort_ind(i)};
% end
% figure; set(gcf,'Position',[100 100 1300 950],'PaperPositionMode','auto')
% barh(R_all0(:,1),0.5,'FaceColor','r')
% hold on
% barh(R_all0(:,2),0.25,'FaceColor','b')
% % legend({'Correlation with UV in exploratory CCA sample','Correlation with UV in replication sample'},'Location','SouthEast')
% % xlabel('Cross-subject correlation with CCA score (UV)'); 
% % title('behaviour CCA correlations');
% set(gca, 'ytick',1:1:length(name1), 'YTickLabel', name1, 'FontSize',17);
% print(gcf,'./loadings/behaviour/behaviour_loading3','-dpng','-r300');