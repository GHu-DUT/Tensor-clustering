clc
clear
close all

tic
%% Parameters need to be changed
load Demodata.mat
data = data-repmat(mean(data'),size(data,2),1)';
runs = 5;
Comp = 2:5;
Method = 'FastICA';
OutPutdir = ['Result_' Method];
%% ICA and Tensor clustering for single subject
f_tensorial_Cluster_Single_Sub(data,runs,Comp,OutPutdir,Method);
%% plot ICA capability scope
switch Method
    case 'FastICA'
        MaxIteration = 100;
    case 'InfomaxICA'
        MaxIteration = 512;
    otherwise
        disp('Unknow method.');
end
f_plot_ICA_Parameter_tensorialClustering(OutPutdir,MaxIteration,Comp);
%%
toc
