function [U,UTV,V,clusterQuality,timeNeeded]=KMeans_TV(pathToSaveResults, numClusters, replicates, tauTV, indexParam, whichData, dataPath)
%Computes the classical K-Means of the dataset X, which has to be
%arranged in such a way, that the observations are ordered rowwise.
%If Replicates = r > 1 and Start is plus (the default), then the software
%selects r possibly different sets of seeds according to the k-means++
%algorithm.

%INPUT:
%pathToSaveResults: Path to save the results of the clustering evaluation.
%numClusters:       Number of Clusters
%replicates:        Number of replicates of the method. If replicates>0,
%                   then the output is the one of the best performed
%                   replicate.
%tauTV:             Regularization parameter for the TV regularizer, which
%                   is applied as a post processing step.
%indexParam:        Index of the vector of the different sigma values for
%                   the parameter tests
%whichData:         Which kind of data should be used? Should be either
%                   'exampleData' or 'ownData' (see README.md).
%dataPath:          In case of 'ownData', the full path to the *.mat file 
%                   with the dataset and all necessary components is
%                   needed (in case of 'exampleData', just use ''; see
%                   README.md.)

%OUTPUT:
%U:                 Cluster Membership Matrix
%UTV:               Cluster membership matrix after TV post processing
%V:                 Matrix, which has the centroids in its rows
%clusterQuality:    Quality of the Clustering in terms of different quality
%                   measures
%timeNeeded:        Needed time for the algorithm in seconds

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.


    %% Choose Dataset
    fprintf('Loading the dataset...');
%     datasetName='colon'; %expected string:
    data=loadData(whichData, dataPath); %load data
    X=data.X;
%     data.name = datasetName;
    fprintf('Done.\n');
    
    %% Set Hyperparameters
    fprintf('Setting the hyperparameters for K-Means...');
    %Set the main hyper parameters
    param = setMainHyperParam(numClusters);
    
    %TV Regularization
    param.tauTV = tauTV;
    param.TVFuzzyFlag = true;
%     param.TVDenoiser = 'Driggs';
    
    param.replicatesKMeans = replicates;
    
    opts = statset('Display','final');
    
    fprintf('Done.\n')
    
    %% Preprocessing Step For Datamatrix
    X=projectNegativeAndHighValues(X, param.epsLim);
    
    %% Start the main Method
    tStart = tic;
    [indexVector,V] = kmeans(X,numClusters,'Distance','sqeuclidean',...
        'Replicates',param.replicatesKMeans,'Options',opts);
    U=indexVectorToClusterMembership(indexVector);

    %% Apply here the TV regularizer and compute clustering quality of U and UTV
    fprintf('Applying TV regularization as post processing step...\n')
    UTV = TV_PostProcessing(U, param.TVFuzzyFlag, param.tauTV,...
        data);

    %% Compute Clustering Quality
    clusterQuality = clusteringQuality(U, data);
    clusterQualityTV = clusteringQuality(UTV, data);

%     clusterQualityVec(end+1) = clusterQualityTV;

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

    warnings.UFuzzy = checkFuzzyClustering(U);
    warnings.UHard = checkHardClustering(convertToHardClustering(U));

    warnings.UTVFuzzy = checkFuzzyClustering(UTV);
    warnings.UTVHard = checkHardClustering(convertToHardClustering(UTV)); %#ok<STRNU>

    %% Compute U^T * U to get the approximative identity matrix
    UU = U'*U; %#ok<NASGU>
    UUTV = UTV'*UTV; %#ok<NASGU>

    timeNeeded=toc(tStart);

    %% Save Everything
    data = rmfield(data,'X'); %#ok<NASGU>
    clearvars -except indexParam pathToSaveResults timeNeeded fileNameAppend numEqualInARowVec updateChangeVec warnings clusterQualityTV UUTV clusterQualityVec i mfilename param U UTV V W clusterQuality i jj data stoppingCause timeAxis UU costFunValues
    nameFile=strcat(pathToSaveResults, '\', mfilename, '_test', num2str(indexParam));
    save(nameFile);
    fprintf('Finished.\n')
end