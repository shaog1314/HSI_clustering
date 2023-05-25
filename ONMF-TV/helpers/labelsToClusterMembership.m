function CMGroundTruth=labelsToClusterMembership(data)
%Computes the cluster membership matrix of the ground truth labels given
%as an image.

%INPUT:
%data:              Struct with the main information relating to the data
%                   without the data itself.

%OUTPUT:
%CMGroundTruth:     Cluster Membership matrix of the ground Truth

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.
    
    NC=max(max(data.labels)); %Number of classes
    CMGroundTruth=zeros(size(data.indexNonzero(:),1), NC);
    
    columnLabels=imageToColumn(data.labels, data);
    
    if any(columnLabels==0)
        error('Some label entries are zero!')
    end
    
    for k=1:NC
        CMGroundTruth(:,k)=(columnLabels==k);
    end
    
end