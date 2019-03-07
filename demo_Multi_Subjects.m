clc
clear
close all

tic
%% Parameters need to be changed
SourceFile = 'data';
OutPutdir = 'GroupICA_FastICA_Result';
runs = 5;
Comp = 2:5;
Method = 'FastICA';
%% ICA decomposition
f_tensorial_Cluster_Multi_Sub(SourceFile,runs,Comp,OutPutdir,Method);
%% plot parameter
switch Method
    case 'FastICA'
        MaxIteration = 100;
    case 'InfomaxICA'
        MaxIteration = 512;
    otherwise
        disp('Unknow method.');
end
f_plot_ICA_Parameter_tensorialClustering(OutPutdir,MaxIteration,Comp);
%% plot Iq for each subject
NumExComp = input(['Please choose the number of extracted components(' num2str(Comp(1)) '-' num2str(max(Comp)) '):']);
load([OutPutdir filesep 'MO_' num2str(NumExComp) filesep 'Iq_AllSub_Temporal.mat']);
load([OutPutdir filesep 'MO_' num2str(NumExComp) filesep 'Matrix_iq.mat']);
figure;
plot(iq,'-*k','linewidth',2);hold on;
plot(Iq_AllSub_Temporal,'o','linewidth',2);grid on;
xlabel('Component#','fontsize',14);
ylabel('Stability index','fontsize',14);
xlim([0.5 NumExComp+0.5]);ylim([0.5 1.1]);
%%
toc
