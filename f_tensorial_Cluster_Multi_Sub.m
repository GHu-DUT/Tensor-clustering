function f_tensorial_Cluster_Multi_Sub(SourceFile,runs,Comp,OutPutdir,Method)
% PURPOSE
% To run ICA decomposition for groupICA
%
% INPUTS
% SourceFile:   (string) the directory of data to be decomposed for each subject
% runs:         (scalar) the number of runs for each model order
% Comp:         (vector) the model order to be selected    
% OutPutdir:    (string) the directory for results to be saved
% Method:       (string) the ICA algorithm to be used: 'FastICA'/'InfomaxICA'

% ver 1.0 030519 GQ

%% Organize data
file = dir([SourceFile filesep '*.mat']);
Data_AllSub_R95 = [];
for isfile = 1:length(file)
    fileName = file(isfile).name;
    disp(['Loading file name:' fileName]);
    struc = load([SourceFile filesep fileName]);
    strucname = fieldnames(struc);
    data = getfield(struc,strucname{1});
    if size(data,1)<size(data,2)
        data = data';
    end
    data = data-repmat(mean(data'),size(data,2),1)';
    [coeff,score,latent] = pca(data);
    Ratio = 0;
    Cont(isfile) = 0;
    while Ratio<0.95
        Cont(isfile) = Cont(isfile)+1;
        Ratio = sum(latent(1:Cont(isfile)))/sum(latent);
    end
    Coeff_AllSub(:,:,isfile) = coeff;
    Data_AllSub_R95 = [Data_AllSub_R95 score(:,1:Cont(isfile))];
end
mkdir(OutPutdir);
save([OutPutdir filesep 'Data_AllSub_R95'],'Data_AllSub_R95','-v7.3');
save([OutPutdir filesep 'Coeff_AllSub'],'Coeff_AllSub','-v7.3');
save([OutPutdir filesep 'Cont'],'Cont','-v7.3');
%%
for isComp = Comp
    ResultFile = [OutPutdir filesep 'MO_' num2str(isComp)];
    mkdir(ResultFile);
    [coeff,score,latent] = pca(Data_AllSub_R95);
    save([OutPutdir filesep 'PCA.mat'],'coeff','score','latent','-v7.3');
    %% ICA decomposition
    switch Method
        case 'FastICA'
            MaxIteration = 100;
            [sR,step]=icassoEst('both', score(:,1:isComp)',runs, 'lastEig', isComp, 'g','tanh', ...
                'approach', 'symm');
        case 'InfomaxICA'
            MaxIteration = 512;
            [sR step]=icassoEst_infomaxICA('both',score(:,1:isComp)',runs, 'lastEig', isComp, 'g', 'tanh', ...
                'approach', 'symm');
        otherwise
            disp('Unknow method.');
    end
    %%
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
    %% Cluster of component matrix
    sR=icassoCluster(sR,0,[],'strategy','AL','simfcn','abscorr','s2d','sim2dis','L','rdim');
    sR=icassoProjection(sR,'cca','s2d','sqrtsim2dis','epochs',75);
    [iq,A,W,S]=icassoShow(sR,'L',isComp,'colorlimit',[.8 .9]);
    save([ResultFile filesep 'Component_A'],'A','-v7.3');
    save([ResultFile filesep 'Component_W'],'W','-v7.3');
    save([ResultFile filesep 'Component_iq'],'iq','-v7.3');
    save([ResultFile filesep 'Component_S'],'S','-v7.3');
    save([ResultFile filesep 'Component_sR'],'sR','-v7.3');
    save([ResultFile filesep 'Component_step'],'step','-v7.3');
    %% Cluster of coefficient matrix
    for isRun = 1:Contruns
        disp(['Calculate features ' num2str(isRun) ' / ' num2str(runs)])
        W = sR.W{isRun};
        S = W*sR.signal;
        Temporal(:,(isRun-1)*isComp+1:isRun*isComp) = coeff(:,1:isComp)/W;
    end
    Feature = rownorm(Temporal');
    SimA = Feature*Feature';
    [iq,A,W,S,sR] = f_Cluster_Feature(SimA,sR,isComp,S);
    save([ResultFile filesep 'Coefficient_iq'],'iq','-v7.3');
    save([ResultFile filesep 'Coefficient_sR'],'sR','-v7.3');
    clear Temporal
    %% Cluster of matrix
    load([ResultFile filesep 'Component_S.mat']);
    SimS = abs(corrw(icassoGet(sR,'W'),icassoGet(sR,'dewhitemat')));
    Similarity = SimA.*SimS;
    [iq,A,W,S,sR] = f_Cluster_Feature(Similarity,sR,isComp,S);
    save([ResultFile filesep 'Matrix_iq'],'iq','-v7.3');
    save([ResultFile filesep 'Matrix_sR'],'sR','-v7.3');
    %% Calculate features
    for isSub = 1:length(file)
        SubsR = sR;
        for isRun = 1:Contruns
            disp(['Subject#' num2str(isSub) ' Calculate features ' num2str(isRun) ' / ' num2str(runs)])
            W = sR.W{1,isRun};
            S = W*sR.signal;
            A = sR.A{isRun};
            SubTemporal = squeeze(Coeff_AllSub(:,1:Cont(isSub),isSub))*...
                coeff(sum(Cont(1:isSub-1))+1:sum(Cont(1:isSub)),1:isComp)/W;
            Temporal(:,(isRun-1)*isComp+1:isRun*isComp) = SubTemporal;  
            SubsR.A(isRun) = {SubTemporal};
        end
        Feature = rownorm(Temporal');
        SimA = Feature*Feature';
        Similarity = SimA.*SimS;
        %% Cluster of matrix
        load([ResultFile filesep 'Component_S']);
        [iq,A,W,S,SubsR] = f_Cluster_Feature(Similarity,SubsR,isComp,S);
        save([ResultFile filesep 'Matrix_Temporal_Iq_Sub#' num2str(isSub)],'iq','-v7.3');
        save([ResultFile filesep 'Matrix_Temporal_sR_Sub#' num2str(isSub)],'SubsR','-v7.3');
        Iq_AllSub_Temporal(:,isSub) = iq;
    end
    save([ResultFile filesep 'Iq_AllSub_Temporal'],'Iq_AllSub_Temporal','-v7.3');
    clear Temporal;
    clear Iq_AllSub_Temporal;
end