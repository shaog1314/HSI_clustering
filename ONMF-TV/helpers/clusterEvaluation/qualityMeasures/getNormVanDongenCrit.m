function VDN=getNormVanDongenCrit(contingencyMatrix)
%Computes the normalized van Dogen criterion (VDN) based on the given
%contingency matrix. 0 is the best value here.

%INPUT:
%contingencyMatrix:     Contingency Matrix 

%OUTPUT:
%VDN:                   Normalized van Dogen criterion (VDN) according to 
%                       Xiong, H.; Li, Z. Clustering Validation Measures.
%                       Data Clustering. Chapman and Hall/CRC, 2014, pp.
%                       571–605.

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.


    n=sum(sum(contingencyMatrix));
    numerator=2*n - sum( max(contingencyMatrix,[],2) ) -...
        sum( max(contingencyMatrix,[],1) );
    denominator=2*n - max(sum(contingencyMatrix,2)) -...
        max(sum(contingencyMatrix,1));
    VDN=numerator/denominator;
    
end