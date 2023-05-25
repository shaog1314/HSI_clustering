function [image]=columnToImage(column, data)
%Computes the image of a given column of a matrix by adding the zero
%spectra.

%INPUT:
%column:        Column of the matrix, which shall be converted into an
%               image
%data:          Struct with the main information relating to the data
%               without the data itself.

%OUTPUT:
%image:         Image of the matrix

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    column=column(:); %Ensure that it is a column
    
    image=zeros(size(data.indexGrid, 1)*size(data.indexGrid, 2), 1);
    image(data.indexNonzero)=column;
    image=reshape(image, size(data.indexGrid));
end