function [column]=imageToColumn(image, data)
%Computes the column of a given image by removing the zero spectra.

%INPUT:
%image:         Image of the matrix
%data:          Struct with all information relating to the data.

%OUTPUT:
%column:        Column of the matrix, which shall be converted into an
%               image

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    column=image(:);
    column=column(data.indexNonzero);
end