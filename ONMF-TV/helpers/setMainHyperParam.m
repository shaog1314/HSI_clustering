function param=setMainHyperParam(numClusters)
%Sets the main hyperparameters for the considered clustering methods. This
%includes:
%   param.iniType
%   param.numClusters
%   param.stoppingCrit
%   param.numEqualInARowMax
%   param.maxIt
%   param.stopLimUV
%   param.rescaling
%   param.epsStabInit
%   param.epsLim

%INPUT:
%numClusters:       Number of clusters, which has to be defined in prior.

%OUTPUT:
%param:             Struct of all hyperparameters of the considered method

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    param.iniType           = 'svdbased'; %svdbased, kmeansPP
    param.numClusters       = numClusters;  %Number of Clusters
    param.stoppingCrit      = 'maxIt'; %Either 'relChange' or 'clusterAssign' or 'PanNg'
    param.numEqualInARowMax = 10; %Maximum number of times that the cluster membership
%                                  assignments do not change, until the flag is set to
%                                  true
    param.maxIt             = 1000;   %Maximal Iteration Number
    param.stopLimUV         = 1e-5;%Stopping Criterion (Relative change) for the matrices U, V.
    param.stopLimUVW        = 1e-5;%Stopping Criterion (Relative change) for the matrices U, V and W.
    param.rescaling         = false; % Flag for Rescaling. Default: false
    param.epsStabInit       = [1e-16, 10^35]; %Constant, which is added at the initialization of the NMF; Before: 1e-8
    param.epsLim            = [1e-16, 10^35]; %Lower and upper projection value during NMF Iterations; Before: 1e-8
    
end