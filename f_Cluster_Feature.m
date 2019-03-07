function [iq,A,W,S,sR] = f_Cluster_Feature(Similarity,sR,Comp,Component_S)
% PURPOSE
% Tensor clustering based multi-information considered similarity matrix
%
% INPUTS
% Similarity:  (matrix) the similarity matrix contain multi-domain information
% sR:          (struct) the ICA structure from ICASSO 
% Comp:        (vector) the model orders    
% Component_S: (matrix) ICA components

% ver 1.0 030519 GQ

sR=icassoCluster(sR,1,Similarity,'strategy','AL','simfcn','abscorr','s2d','sim2dis','L','rdim');
sR=icassoProjection(sR,'cca','s2d','sqrtsim2dis','epochs',75);
[iq,A,W,S]=icassoShow(sR,'L',Comp,'colorlimit',[.8 .9]);
rho = corr(S',Component_S');
[~,order] = max(abs(rho)');
S = S(order,:);
iq = iq(order);
W = W(order,:);
A = A(:,order);
sR.order = order;

end