clear all; close all; clc

% Based on https://git.fmrib.ox.ac.uk/falmagro/ukb_unconfound_v2/-/blob/master/conf_processing/common_matlab/duplicateDemedianNormBySite.m
% and https://git.fmrib.ox.ac.uk/falmagro/ukb_unconfound_v2/-/blob/master/conf_processing/common_matlab/duplicateCategorical.m

addpath('/mnt/data0/home/khu/CCA/MatlabScripts/FSLNets')
addpath('/mnt/data0/home/khu/CCA/MatlabScripts')
INPUT = 'eid_train.csv';
% INPUT = 'eid_test.csv';
% INPUT = 'eid_19157.csv';

%% Load confound information
S = load(sprintf('%s/%s','/mnt/data0/home/khu/CCA/Data/mental_sMRI/eid_19157/random4',INPUT));
% S = importdata(sprintf('%s/%s','/mnt/data0/home/khu/CCA/Data/mental_sMRI/eid_19157',INPUT));
% S = S.data;
S(8608) = [];
S(10444) = [];
Conf = readtable(['/mnt/data0/home/khu/CCA/Data/mental_sMRI/confound_19175.csv'],'FileType','text');
[~,s,~] = intersect(table2array(Conf(:,1)),S); Conf = Conf(s,:);
H = get_UKB_headers(Conf);

%% Create site varriables
finalConfs = [];
n1 = strfind(H,'54-2.0'); n1 = find(~cellfun(@isempty,n1));
site = table2array(Conf(:,n1)); clear n1
values = site;
diffValues = unique(values);
numDiffValues = length(diffValues);
if numDiffValues > 1
    numNewCols = numDiffValues-1;
    newConfound = zeros(size(Conf,1), numDiffValues-1);
    
    % Subjects (of this site) with the first value get a -1 in all new columns
    indFinal = find(values == diffValues(1));
    newConfound(indFinal,:)=-1;
    
    % Each new value gets a new column. Subjects (of this Site) with
    % this value have 1 in this column. All other subjects have a 0.
    for j=2:numDiffValues
        indFinal = find(values == diffValues(j));
        newConfound(indFinal,j-1)=1;
        
        indNotZero = find(newConfound(:,j-1) ~= 0);
        tmpVar = newConfound(indNotZero,j-1);
        
        % Normalising only makes sense if there is more than 1
        % different value
        newConfound(indNotZero,j-1) = nets_normalise(newConfound(indNotZero,j-1));
    end
    
end

%% Creating other variables
% Age
n1 = strfind(H,'21003-2.0'); n1 = find(~cellfun(@isempty,n1));
values = table2array(Conf(:,n1)); clear n1
% Age squared
values(:,2) = values(:,1).^2;
% Sex
n1 = strfind(H,'31-0.0'); n1 = find(~cellfun(@isempty,n1));
values = [values table2array(Conf(:,n1))]; clear n1
% Age * sex
values(:,4) = values(:,1).*values(:,3);
% Head size
n1 = strfind(H,'25000-2.0'); n1 = find(~cellfun(@isempty,n1));
values = [values table2array(Conf(:,n1))]; clear n1
% Date of attending assessment centre
n1 = strfind(H,'53-2.0'); n1 = find(~cellfun(@isempty,n1));
days = table2array(Conf(:,n1)); clear n1
days = datenum(days); days = days - (min(days)-1);
values = [values days]; clear days
% Date squared
values(:,7) = values(:,6).^2;
% interval years (Date of attending assessment centre - Date of completing mental health questionnaire)
n1 = strfind(H,'53-2.0'); n1 = find(~cellfun(@isempty,n1));
days1 = table2array(Conf(:,n1)); 
days1 = datenum(days1); 
n2 = strfind(H,'20400-0.0'); n2 = find(~cellfun(@isempty,n2));
days2 = table2array(Conf(:,n2)); 
days2 = datenum(days2); 
interval = days1 - days2;
values = [values interval]; clear days1; clear days2; clear n1; clear n2;


%% Split regressors for sites
numVars  = size(values,2);
subjectIndicesBySite{1} = find(site==11025);
subjectIndicesBySite{2} = find(site==11026);
subjectIndicesBySite{3} = find(site==11027);
% subjectIndicesBySite{4} = find(site==11028);
numSites = length(subjectIndicesBySite);
% sum(site==11028)

% Demedian globally each column
for i = 1 : numVars
    values(:,i) = values(:,i) - nanmedian(values(:,i));
    madV = mad(values(:,i),1) * 1.48;
    if madV < eps
        madV = nanstd(values(:,i));
    end
    values(:,i) = values(:,i) / madV;
end

% Split the colums: 1 per site. 0-padding by default.
finalConfs = [finalConfs zeros(size(Conf,1),numSites*numVars)];

% For each Site
for i = 1:numSites
    % For each variable (column) that is received in "values"
    for j = 1:numVars
        V = values(subjectIndicesBySite{i},j);
        
        medTP = nanmedian(V);
        
        V(find(V > 8)) =  NaN;
        V(find(V <-8)) =  NaN;
        
        V(isnan(V)) = medTP;
        if any(V) == 1 & length(V) > 1% ���ڷ�0Ԫ�أ���ȫΪ0��any���������������Ƿ��з���Ԫ�أ�����У��򷵻�?1�����򣬷���0
            V = nets_normalise(V);
        end
        finalConfs(subjectIndicesBySite{i}, (numVars*(i-1) + j)) = V;
    end
end

conf=[finalConfs, newConfound];
save(sprintf('./eid_19157/random4/new_15323/Confounds_Subjects_%s.mat',INPUT(5:end-4)), 'conf');
dlmwrite(sprintf('./eid_19157/random4/new_15323/Confounds_Subjects_%s.csv',INPUT(5:end-4)),conf,'precision',25);

