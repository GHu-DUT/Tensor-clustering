function f_plot_ICA_Parameter_tensorialClustering(ResultName,MaxIteration,Comp)
% PURPOSE
% To plot ICA parameters based Tensor clustering
%
% INPUTS
% ResultName:   (string) results saved directory
% MaxIteration: (scalar) max iteration defined for algorithm
% Comp:         (vector) the range of number of extracted components
%

% ver 1.0 073019 GQ

maxRun = 0;
for isComp = Comp
    isComp
    load([ResultName filesep 'MO_' num2str(isComp) filesep 'Matrix_iq.mat']);
    tmp = load([ResultName filesep 'MO_' num2str(isComp) filesep 'Component_step.mat']);
    tmp = tmp.step;
    load([ResultName filesep 'MO_' num2str(isComp) filesep 'Matrix_sR.mat']);
    runs = length(sR.W);
    maxRun = max([maxRun runs]);
    partition = sR.cluster.partition;
    [Iq, A, W, S, index2centrotypes]=icassoResult(sR,isComp);
    for isRun = 1:runs
        Order = partition(isComp,(isRun-1)*isComp+1:isRun*isComp);
        W_isRun = sR.W{isRun};
        W_isRun(Order,:) = W_isRun;
        rho_W(isRun) = mean(abs(diag(corr(W_isRun',W'))));
        A_isRun = sR.A{isRun};
        A_isRun(:,Order) = A_isRun;
        rho_A(isRun) = mean(abs(diag(corr(A_isRun,A))));
    end
    rho(isComp) = mean(rho_A);
    Patameter(isComp,1) = mean(iq);
    Patameter(isComp,2) = std(iq);
    Patameter(isComp,3) = mean(tmp(tmp<MaxIteration));
    Patameter(isComp,4) = std(tmp(tmp<MaxIteration));
    Patameter(isComp,5) = size(tmp(tmp<MaxIteration),2);
end
%%
maxComp = max(Comp);
figure;
set(gca,'fontsize',14);
plot(rho,'ok','linewidth',2);grid on;
ylabel('Similairty of coefficient matrix');
xlabel('Number of extracted components');
figure
plot(Patameter(1:maxComp,3),'xg','linewidth',2);hold on
plot(Patameter(1:maxComp,4),'or','linewidth',2);hold off
set(gca,'fontsize',14);
ylabel('Number of steps')
xlabel('Number of extracted components');
ylim([-10 MaxIteration+40]);
set(gca,'fontsize',14);
legend(['Mean of numbers of used steps' sprintf('\n') 'to converge among ' int2str(runs) ' runs'],...
    ['SD of numbers of used steps' sprintf('\n') 'to converge among ' int2str(runs) ' runs'],'Location','Best');
grid on;xlim([min(Comp)-0.5 max(Comp)+1]);
figure;
plot(Patameter(1:maxComp,1),'+m','linewidth',2);hold on
plot(Patameter(1:maxComp,2),'*r','linewidth',2);hold off
ylim([-0.1 1.1]);
set(gca,'fontsize',14);
xlabel('Number of extracted components');
ylabel('Magnitude of Iq from tensor clustering')
set(gca,'fontsize',14);
legend('Mean of Iqs of extracted components','SD of Iqs of extracted components','Location','Best');
grid on;xlim([min(Comp)-0.5 max(Comp)+1]);
figure;
plot(Patameter(1:maxComp,5),'v','linewidth',2);
ylim([-5 maxRun+5]);xlim([min(Comp)-0.5 max(Comp)+1]);
set(gca,'fontsize',14);
xlabel('Number of extracted components');
set(gca,'fontsize',14);
ylabel('Number of convergenced runs');
grid on
end