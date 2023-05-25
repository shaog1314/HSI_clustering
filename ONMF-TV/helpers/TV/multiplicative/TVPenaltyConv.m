function [TVU, gradU3D]=TVPenaltyConv(U, data, param)
%Calculates the smoothed, discrete isotropic TV penalty term of a given
%matrix U, which has images in its columns without the zero spectra.
%The weights of the TV penalty term (Psi) is assumed here to be one for all
%entries.

%INPUT
%U:                 Cluster Membership Matrix U
%data:              Struct with all information relating to the data.
%param:             struct with all hyperparameters

%OUTPUT
%TVU:               TV Penalty of the matrix U
%gradU3D:           Gradient of all images of the columns of U in a 3-way
%                   array.

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.


if nargin ~= 3
    error('Wrong number of input arguments.')
end


gradientU = zeros(size(U));
gradU3D = zeros(size(data.indexGrid,1), size(data.indexGrid,2), size(U,2));

%Define the Kernels for the convolutions
deltaXKernel=[-1; 1];
deltaYKernel=[-1 1];

for i=1:size(U,2)
    
    UImage = zeros(size(data.indexGrid,1)*size(data.indexGrid,2),1);
    UImage(data.indexNonzero) = U(:,i);
    UImage = reshape(UImage, size(data.indexGrid));
    
%     gradientUImage = zeros(size(data.indexGrid));
    
    %deltaX
    deltaX = conv2(UImage,deltaXKernel,'same');
    deltaX(end,:) = 0;

    %deltaY
    deltaY = conv2(UImage,deltaYKernel,'same');
    deltaY(:,end) = 0;
    
    gradientUImage = sqrt(deltaX.*deltaX + deltaY.*deltaY +...
        param.epsilonTV^2*ones(size(data.indexGrid)));
    
    gradU3D(:,:,i) = gradientUImage;
    
    %Reshape the gradient of the image to a column
    gradientUImage = imageToColumn(gradientUImage, data);
    
%     gradientUImage=reshape(gradientUImage, [size(data.indexNonzero,1)*...
%         size(data.indexNonzero,2), 1]);
%     
%     gradientUImage=gradientUImage(data.indexNonzero);
    
    gradientU(:,i) = gradientUImage;
    
end

    
%Final calculation of the TV penalty
    Psi=repmat(param.TVWeights, size(U,1), 1);
    TVU = sum(sum(Psi.*gradientU));