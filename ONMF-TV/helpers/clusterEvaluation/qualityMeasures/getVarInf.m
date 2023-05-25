function VI=getVarInf(contingencyMatrix)
%Computes the variational of Information (VI) based on the given
%contingency Matrix. 0 is the best value.

%INPUT:
%contingencyMatrix:     Contingency Matrix 

%OUTPUT:
%VI:                    Variational of information (VI)
%                       according to 
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
    PiPj=sum(P,2)*sum(P,1);%This is a matrix of the same size as P.
    
    logMatrix1=sum(P,2).*log2(sum(P,2));
    logMatrix1(isnan(logMatrix1))=0; %Take care here of the NaN entries
    logMatrix2=sum(P,1).*log2(sum(P,1));
    logMatrix2(isnan(logMatrix2))=0; %Take care here of the NaN entries
    logMatrix3=P.*log2( P./PiPj );
    logMatrix3(isnan(logMatrix3))=0; %Take care here of the NaN entries
    
    VI=-sum( logMatrix1 ) - sum( logMatrix2 )-2*sum(sum( logMatrix3 ));
    
end