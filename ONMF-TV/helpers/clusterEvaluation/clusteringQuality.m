function clusterQuality=clusteringQuality(U, data)
%Compute all the considered quality measures for the given hard clustering
%result given by the cluster membership matrix 'clustering'. See also Data 
%Clustering - Algorithms and Applications: Clustering Validation Measures, 
%Hui, Zhongmou.

%INPUT:
%U:                 Cluster membership matrix with the result of the used
%                   algorithm
%data:              Struct with all information relating to the data. 

%OUTPUT:
%clusterQuality:    Struct, which contains every considered clustering
%                   quality measure for the given hard clustering result 
%                   given by the cluster membership matrix 'clustering'.

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    clustering = convertToHardClustering(U);

    contingencyMatrix=getContingencyMatrix(clustering, data);
    
    clusterQuality.E=getEntropy(contingencyMatrix);
    clusterQuality.P=getPurity(contingencyMatrix);
    clusterQuality.VI=getVarInf(contingencyMatrix);
    clusterQuality.VD=getVanDongenCrit(contingencyMatrix);
    clusterQuality.VIN=getNormVarInf(contingencyMatrix);
    clusterQuality.VDN=getNormVanDongenCrit(contingencyMatrix);
    
    %Get the distance of U to the identity matrix
    [clusterQuality.distUtoIl1, clusterQuality.distUtoIl2,...
        clusterQuality.distUtoImax]=getDistanceToIdentity(U);
    
end