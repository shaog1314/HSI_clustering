function [gradImage]=gradientImage(image, epsTV)
%Computes the gradient of a given image based on the approximative
%differentiable, isotropic Total Variation penalty term with Neumann
%boundary conditions. (See for that L. Condat, Discrete Total Variation:
%New Definition and Minimization, 2017.

%Careful: The images in U have to include the (possible) zero-spectra.

%INPUT:
%image:             Input image without Zero-Spectra(!)
%epsTV:             Small constant to ensure the differentiability of the
%                   TV penalty term

%OUTPUT:
%gradImage:         Gradient of the Total Variation penalty term.

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

%Define the kernels:
deltaXKernel=[1; 0; -1]/2;
deltaYKernel=[1 0 -1]/2;
deltaXXKernel=[1; -2; 1];
deltaYYKernel=[1 -2 1];
deltaXYKernel=[1 0 -1; 0 0 0; -1 0 1]/4;

%Compute the different convolutions, i.e. the discretizations of the
%derivatives and apply the Neumann boundary conditions.

%deltaX
deltaX = conv2(image,deltaXKernel,'same');
deltaX([1; end],:) = 0;

%deltaY
deltaY = conv2(image,deltaYKernel,'same');
deltaY(:,[1,end]) = 0;

%deltaXX
deltaXX = conv2(image,deltaXXKernel,'same');
deltaXX([1; end],:) = 0;

%deltaY
deltaYY = conv2(image,deltaYYKernel,'same');
deltaYY(:,[1,end]) = 0;

%deltaXY
deltaXY = conv2(image,deltaXYKernel,'same');
deltaXY([1; end],:) = 0;
deltaXY(:,[1,end]) = 0;

%Compute gradTV based on the computed discretization
gradImage = (epsTV*(deltaXX + deltaYY) + ((deltaX.*deltaX).*deltaYY) + ...
    ((deltaY.*deltaY).*deltaXX) - 2.*deltaX.*deltaY.*deltaXY ) ./ (sqrt(...
    (deltaX.*deltaX) + (deltaY.*deltaY) + epsTV*ones(size(image)) )).^3;

end