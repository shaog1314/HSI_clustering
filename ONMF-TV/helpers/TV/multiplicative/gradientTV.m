function [gradTV]=gradientTV(U, data, epsTV)
%Gathers together the gradient of the images computed by gradientImage.m
%based on the approximative differentiable, isotropic Total Variation
%penalty term with Neumann boundary conditions. (See for that L. Condat,
%Discrete Total Variation: New Definition and Minimization, 2017.)

%Careful: The images in U have to be without the (possible) zero-spectra.

%INPUT:
%U:                 Cluster Membership Matrix without Zero-Spectra(!)
%data:              Struct with all information relating to the data. 
%epsTV:             Small constant to ensure the differentiability of the
%                   TV penalty term

%OUTPUT:
%gradTV:            Gradient of the Total Variation penalty term for the
%                   given matrix U, which has the same size as U.

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

gradTV = zeros(size(U));

for k=1:size(U,2)
    image = columnToImage(U(:,k), data);
    image = reshape(image, size(data.indexGrid));
    gradImage=gradientImage(image, epsTV);
    gradImage=gradImage(:);
    gradImage=gradImage(data.indexNonzero);
    gradTV(:,k)=gradImage;
end


end