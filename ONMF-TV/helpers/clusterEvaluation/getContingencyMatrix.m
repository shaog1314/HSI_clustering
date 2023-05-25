function contMatrix=getContingencyMatrix(CM, data)
%Computes the contingency matrix based on the computed hard clustering 
%result and the ground Truth of the dataset. Columnwise the true class
%labels and row wise the computed partitions.

%INPUT:
%CM:                Cluster membership matrix with the hard clustering 
%                   result of the used algorithm
%data:              Struct with the main information relating to the data
%                   without the data itself. 

%OUTPUT:
%contMatrix:        Contingency Matrix

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.
    
    if nargin~=2
        error('Wrong number of input arguments')
    end
    %Convert possible fuzzy clustering to hard clustering
    CM=convertToHardClustering(CM);
    
    %Set all the nonzero entries to one. This is important later on to
    %count the entities in the respective partitions.
    CM(CM>0)=1;
%     contMatrix=zeros(size(clustering,1),max(labels(:)));
    
    CMGroundTruth=labelsToClusterMembership(data);
    
    if size(CM)~=size(CMGroundTruth)
        error('The sizes of the cluster membership matrices of the ground Truth and the clustering result are not equal');
    end
    
    contMatrix=CM'*CMGroundTruth;
    
    %Test if sum of contingency matrix is equal to the number of nonzero
    %elements in the labels
    if sum(sum(contMatrix))~=sum(sum(data.labels>0))
        error('Sum of contingency matrix is not equal to the number of nonzero elements in the labels');
    end
    
end