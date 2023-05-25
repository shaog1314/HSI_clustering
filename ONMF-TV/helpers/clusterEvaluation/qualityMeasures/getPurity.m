function purity=getPurity(contingencyMatrix)
%Computes the purity based on the given contingency Matrix. 1 is the best
%value.

%INPUT:
%contingencyMatrix:     Contingency Matrix 

%OUTPUT:
%purity:                Purity according to 
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

    P=contingencyMatrix/sum(sum(contingencyMatrix));
    Pij_Pi=P./( repmat( sum(P,2),1,size(P,2) ) ); %Matrix of the same size as P
%     PiPj=sum(P,2)*sum(P,1);%This is a matrix of the same size as P.
    purity=sum( sum(P,2) .* ( max(Pij_Pi,[],2) ) );
    
end