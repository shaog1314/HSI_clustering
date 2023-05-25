function [indexNonzero] = getIndexNonzeroFromIndexGrid(indexGrid)
%Get the indexNonzero vector based on the given indexGrid of the dataset.

%INPUT:
%indexGrid:         The indexGrid as a matrix of the dataset

%OUTPUT:
%indexNonzero:      The vector containing the index values of the non-zero
%                   spectra

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    %Convert to double precision
    indexGrid=double(indexGrid);
    
    indexNonzero = indexGrid(:);
    indexNonzero(indexNonzero==0)=[];
    
    t=0;
    for j=1:size(indexGrid, 2)
        for i=1:size(indexGrid, 1)
            if indexGrid(i, j)>0
                t=t+1;
                indexNonzero(t) = (j-1)*size(indexGrid,1)+i;
            end
        end
    end

end