function indexVector=clusterMembershipToIndexVector(U)
%Computes the indexVector based on the (fuzzy) cluster membership matrix U,
%which is given without the zero spectra. The indexVector also does not
%contain any zero spectra.

%INPUT:
%U:                 (Fuzzy) Cluster membership matrix (without the
%                   zero-spectra)

%OUTPUT:
%indexVector:       indexVector of the hard clustering, which was given by
%                   the cluster membership matrix (without the zero 
%                   spectra)

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    Uhard=convertToHardClustering(U);
    
    indexVector = zeros(size(Uhard,1),1);
    
    for i=1:size(Uhard,2)
        indexVector = indexVector + Uhard(:,i)*i;
    end
    
    if any(indexVector==0)
        error('Some index entries are zero!')
    end
    
end