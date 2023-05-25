function CM=indexVectorToClusterMembership(indexVector)
%Computes the cluster membership matrix of the index vector obtained by a
%clustering algorithm. The index vector should not contain the zero
%spectra.

%INPUT:
%indexVector:       Index Vector of the hard clustering

%OUTPUT:
%CM:                Cluster Membership matrix of the hard clustering

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.
    
    %Ensure, that the indexVector is a column vector
    indexVector=indexVector(:);

    NP=max(indexVector); %Number of partitions
    CM=zeros(size(indexVector(:),1), NP);
    
    if any(indexVector==0)
        error('Some index entries are zero!')
    end
    
    for k=1:NP
        CM(:,k)=(indexVector==k);
    end
    
end