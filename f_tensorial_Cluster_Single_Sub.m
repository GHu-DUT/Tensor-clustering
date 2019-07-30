function f_tensorial_Cluster_Single_Sub(data,runs,Comp,OutPutdir,Method)
% PURPOSE
% To run ICA decomposition for single subject
%
% INPUTS
% data:       (matrix) the data to be decomposed
% runs:       (scalar) the number of runs for each model order
% Comp:       (vector) the model order to be selected    
% OutPutdir:  (string) the directory for results to be saved
% Method:     (string) the ICA algorithm to be used: 'FastICA' / 'InfomaxICA'

% ver 1.0 073019 GQ


for isComp = Comp
    ResultFile = [OutPutdir filesep 'MO_' num2str(isComp)];
    mkdir(ResultFile);
    %% PCA
    [coeff,score,latent] = pca(data);
    save([OutPutdir filesep 'PCA.mat'],'coeff','score','latent');
    %% ICA decomposition
    switch Method
        case 'FastICA'
            MaxIteration = 100;
             [sR,step]=icassoEst('both',score(:,1:isComp)',runs, 'lastEig', isComp, 'g','tanh', ...
                 'approach', 'symm');
        case 'InfomaxICA'
            MaxIteration = 512;
            [sR step]=icassoEst_infomaxICA('both',score(:,1:isComp)',runs, 'lastEig', isComp, 'g', 'tanh', ...
                'approach', 'symm');
        otherwise
            disp('Unknow method.');
    end
    %% Reject runs don't converged
    Contruns = sum(step<MaxIteration);
    New_sR.mode = sR.mode;
    New_sR.signal = sR.signal;
    New_sR.fasticaoptions = sR.fasticaoptions;
    New_sR.index = sR.index(1:Contruns*isComp,:);
    New_sR.A = sR.A(step<MaxIteration);
    New_sR.W = sR.W(step<MaxIteration);
    New_sR.whiteningMatrix = sR.whiteningMatrix;
    New_sR.dewhiteningMatrix = sR.dewhiteningMatrix;
    New_sR.cluster = sR.cluster;
    New_sR.projection = sR.projection;
    sR = New_sR;
    %% Cluster for component
    sR=icassoCluster(sR,0,[],'strategy','AL','simfcn','abscorr','s2d','sim2dis','L','rdim');
    sR=icassoProjection(sR,'cca','s2d','sqrtsim2dis','epochs',75);
    [iq,A,W,S]=icassoShow(sR,'L',isComp,'colorlimit',[.8 .9]);
    save([ResultFile filesep 'Component_step'],'step','-v7.3');
    save([ResultFile filesep 'Component_A'],'A','-v7.3');
    save([ResultFile filesep 'Component_W'],'W','-v7.3');
    save([ResultFile filesep 'Component_iq'],'iq','-v7.3');
    save([ResultFile filesep 'Component_S'],'S','-v7.3');
    save([ResultFile filesep 'Component_sR'],'sR','-v7.3');
    %% Cluster of coefficient
    for isRun = 1:Contruns
        disp(['Calculate features ' num2str(isRun) ' / ' num2str(runs)])
        A = sR.A{isRun};
        Temporal(:,(isRun-1)*isComp+1:isRun*isComp) = A;
    end
    Feature = rownorm(Temporal');
    SimA = Feature*Feature';
    [iq,A,W,S,sR] = f_Cluster_Feature(SimA,sR,isComp,S);
    save([ResultFile filesep 'Coefficient_iq'],'iq','-v7.3');
    save([ResultFile filesep 'Coefficient_sR'],'sR','-v7.3');
    %% Tensorial clustering
    SimS = abs(corrw(icassoGet(sR,'W'),icassoGet(sR,'dewhitemat')));
    Similarity = SimA.*SimS;
    load([ResultFile filesep 'Component_S.mat']);
    [iq,A,W,S,sR] = f_Cluster_Feature(Similarity,sR,isComp,S);
    save([ResultFile filesep 'Matrix_iq'],'iq','-v7.3');
    save([ResultFile filesep 'Matrix_sR'],'sR','-v7.3');
    clear Temporal
end