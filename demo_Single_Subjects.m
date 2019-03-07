clc
clear
close all

tic
%% Parameters need to be changed
load megdata.mat
megdata = megdata-repmat(mean(megdata'),size(megdata,2),1)';
runs = 5;
Comp = 2:10;
Method = 'FastICA';
OutPutdir = ['Result_' Method];
%% ICA and Tensor clustering for single subject
f_tensorial_Cluster_Single_Sub(megdata',runs,Comp,OutPutdir,Method);
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
