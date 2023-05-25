function [UTV,clusterQuality,clusterQualityTV,clusterQualityVecNew,warnings,UU,UUTV]=tasksAfterMainAlgForSepMethods(stoppingCause,param,U,warningsOld,data,clusterQualityVec)
%Performs all tasks (except the saving) after the main algorithm of the
%separated methods, which includes the TV regularization step, which is
%performed as a post processing step.

%INPUT:
%stoppingCause:         Cause of the stopping of the main algorithm. Can be
%                       'relChange', clusterAssign' or 'maxIt'
%param:                 Struct of all hyperparameters of the considered
%                       method
%U:                     Cluster membership matrix
%warningsOld:           Struct, which consists of the warnings of the
%                       initialization of U
%data:                  Struct with all information relating to the data.
%clusterQualityVec:     multidimensional struct, which contains the cluster
%                       quality measures from all iterations steps.

%OUTPUT:
%UTV:                   Cluster membership matrix after TV post processing
%clusterQuality:        Struct, which contains the cluster quality measures
%                       of U
%clusterQualityTV:      Struct, which contains the cluster quality measures
%                       of UTV
%clusterQualityVec:     Updated clusterQualityVec
%warnings:              Struct, which consists of all possible warnings
%                       which can occur during the considered method.
%UU:                    Product U'*U
%UUTV:                  Product UTV'*UTV


% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

%% Output of the stopping cause
    
    switch stoppingCause
        case 'relChange'
            fprintf('Stopping Cause: Relative Change of the matrices U, V and W.\n')
        case 'clusterAssign'
            fprintf('Stopping Cause: No change of cluster assignments during %d iterations.\n',...
                param.numEqualInARowMax)
        case 'maxIt'
            fprintf('Stopping Cause: Maximum number of iterations (%d) reached.\n',param.maxIt)
    end
    
    %% Apply here the TV regularizer and compute clustering quality of U and UTV
    fprintf('Applying TV regularization as post processing step...\n')
    UTV = TV_PostProcessing(U, param.TVFuzzyFlag, param.tauTV,...
        data);
    
    %% Compute Clustering Quality
    clusterQuality = clusteringQuality(U, data);
    clusterQualityTV = clusteringQuality(UTV, data);
    
    clusterQualityVecNew = clusterQualityVec;
    clusterQualityVecNew(end+1) = clusterQualityTV;
    
    msg=['Cluster Quality after the TV post processing:',...
        num2str(clusterQualityTV.distUtoIl1), '\nDistance UTV to identity (l2): ',...
          num2str(clusterQualityTV.distUtoIl2), '\nDistance UTV to identity (max. Value): ',...
          num2str(clusterQualityTV.distUtoImax), '\n\nPurity: ',...
          num2str(clusterQualityTV.P), '\nEntropy: ', num2str(clusterQualityTV.E),...
          '\nVIN: ', num2str(clusterQualityTV.VIN), '\nVDN: ',...
          num2str(clusterQualityTV.VDN),'\n'];
    fprintf(msg)
    
    fprintf('Post processing with TV done! Last tasks and saving procedure...\n')
    
    %% Check Clusterings U and UTV
    warnings = warningsOld;
    
    warnings.UFuzzy = checkFuzzyClustering(U);
    warnings.UHard = checkHardClustering(convertToHardClustering(U));
    
    warnings.UTVFuzzy = checkFuzzyClustering(UTV);
    warnings.UTVHard = checkHardClustering(convertToHardClustering(UTV));
    
    %% Compute U^T * U to get the approximative identity matrix
    UU = U'*U;
    UUTV = UTV'*UTV;
    
end