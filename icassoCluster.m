function sR=icassoCluster(sR,ClusterMode,Similarity,varargin)
% demo sR=icassoCluster(sR,1,Similarity,'strategy','AL','simfcn','abscorr','s2d','sim2dis','L','rdim');
% Init some variables

% total number of estimates
M=icassoGet(sR,'M');

%reduced data dimension
rdim=icassoGet(sR,'rdim');

% Set default parameters
default={'simfcn','abscorr','s2d','sim2dis','strategy','AL','L','rdim'};

%% Check optional arguments and add defaults
clusterparameters=processvarargin(varargin,default);
num_of_args=length(clusterparameters);

%% check arguments
for i=1:2:num_of_args;
    switch lower(clusterparameters{i})
        
        case 'simfcn'
            simfcn=clusterparameters{i+1};
            
            % Explicit similarity matrix?
            if isnumeric(simfcn),
                if size(simfcn,1)==M & size(simfcn,2)==M,
                    sR.cluster.similarity=simfcn;
                    sR.cluster.simfcn='<similarities given explicitly>';
                else
                    error('Explicitly given similarity matrix has wrong size!');
                end
            else
                % should be a string
                switch lower(simfcn)
                    case 'abscorr'
                        ; % ok
                        sR.cluster.simfcn=lower(simfcn);
                    otherwise
                        error('''simfcn'' must be string ''abscorr'' or an MxM similarity matrix');
                end
            end
        case 's2d'
            s2dfcn=lower(clusterparameters{i+1});
            if ~ischar(s2dfcn),
                error('''s2d'' must be a string (name of a function)');
            end
            sR.cluster.s2d=s2dfcn;
        case 'l'
            L=clusterparameters{i+1};
            if isnumeric(L),
                % The user has specified max number for clusters
                
                % Check L
                if fix(L)~=L,
                    error('''L'' must be an integer.');
                elseif L<2,
                    error('''L'' must be at least 2.');
                elseif L>M,
                    error('''L'' cannot be more than the number of estimates.');
                end
            else
                if ~strcmp(lower(L),'rdim'),
                    error('''L'' expects an integer value or ''rdim''.');
                end
                % set (reduced) data dimension
                L=icassoGet(sR,'rdim');
            end
            
            if L>100,
                warning(['R-index requested for more that 100 clusters: this can' ...
                    ' be heavy...']);
            end
            
        case 'strategy'
            strategy=clusterparameters{i+1};
            if ~ischar(strategy),
                error('''strategy'' must be a string');
            end
            
            % we are case insensitive
            strategy=upper(strategy);
            sR.cluster.strategy=strategy;
            
            switch sR.cluster.strategy
                case {'AL','CL','SL'}
                    ; % hierarchical clustering
                otherwise
                    error(['Strategy ' strategy ' not implemented.']);
            end
        otherwise
            error(['Indentifier ' clusterparameters{i} ' not recognized.']);
    end
end

%%%% Compute similarities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(sR.cluster.simfcn)
    case '<similarities given explicitly>'
        ; % already handled
    case 'abscorr'
        if ClusterMode==1
            sR.cluster.similarity = abs(Similarity);
        else
            sR.cluster.similarity=abs(corrw(icassoGet(sR,'W'),icassoGet(sR,'dewhitemat')));
        end
        %just to make sure
        sR.cluster.similarity(sR.cluster.similarity>1)=1;
        sR.cluster.similarity(sR.cluster.similarity<0)=0;
end
%%%%% Convert to dissimilarities using .s2d
D=feval(sR.cluster.s2d, sR.cluster.similarity);
%%%% Make partition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sR.cluster.partition,sR.cluster.dendrogram.Z,sR.cluster.dendrogram.order]=...
    hcluster(D,sR.cluster.strategy);
%%%%% Compute cluster validity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sR.cluster.index.R=ones(M,1)*NaN;
sR.cluster.index.R(1:L,1)=rindex(D,sR.cluster.partition(1:L,:));
