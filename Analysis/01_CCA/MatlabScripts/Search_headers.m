clear all; close all; clc
load UKB_headers_all.mat

strcmp('54-2.0',H41)
find(ans==1)
strncmp('20400',H41,5)
find(ans==1)

for n = 1:length(sarah)
    IN = [num2str(sarah(n,1)) '-'];
    F = strncmp(IN,H41,length(IN));
    sarah(n,2) = sum(F); clear F
end
for n = 1:length(jonathan)
    IN = [num2str(jonathan(n,1)) '-'];
    F = strncmp(IN,H41,length(IN));
    jonathan(n,2) = sum(F); clear F
end