function UPost=TV_PostProcessing(U, TVFuzzyFlag, lambda, data)
%Performs a TV denoising as post processing on the cluster membership
%matrix U.

%INPUT:
%U:                 Cluster membership matrix
%TVFuzzyFlag:       Flag. True, if the TV regularization shall be applied
%                   to the fuzzy cluster membership matrix U. False, if it
%                   is applied to the hard clustering.
%lambda:            Regularization parameter of the TV denoising algorithm
%data:              Struct with all information relating to the data.

%OUTPUT:
%UPost:             Post processed Cluster Membership Matrix

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    
    UPost = 100*ones(size(U));
%     niter = 100; %Or higher.
    
    if ~TVFuzzyFlag
        fprintf('TV Postprocessing on hard U with Driggs...')
        U = convertToHardClustering(U);
    else
        fprintf('TV Postprocessing on fuzzy U with Driggs...')
    end
    for ii=1:size(U,2)
        UPost(:,ii) = proxTV(U(:,ii), data, lambda);
    end
    
    %Apply projection step to avoid negative entries in UTV
    UPost(UPost<0)=0;
end

function [columnUOut]   =   proxTV(columnUIn, data, tauTV)

    param.verbose = 0;
    param.maxit = 100;
    imageU = columnToImage(columnUIn, data);
%     noiseLevel = 1e-4;
%     imageU = max(imageU+max(abs(imageU(:)))*noiseLevel*...
%         randn(size(imageU)),0);
    [imageUTV,~] = prox_tv_19(imageU, tauTV, param);
    columnUOut = imageToColumn(imageUTV, data);
    
end