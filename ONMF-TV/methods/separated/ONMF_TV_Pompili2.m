function [U,UTV,V,clusterQuality,timeNeeded]=ONMF_TV_Pompili2(pathToSaveResults, numClusters, tauTV, alpha0Inp, stepBetaInp, iniType, iniTypePompili, indexParam, whichData, dataPath)
%Computes the (orthogonal) NMF based on an augmented lagrangian method for
%ONMF in Pompili et al., 2014 (Neurocomputing) (ONP-MF). The
%datamatrix has the dimensions MxN, where M denotes the number of
%observations and N the number of features.
%There are hyperparameters which have to be chosen a priori, but we will
%leave it with the default case. We use the proposed stopping criterion of
%Pompili et al., 2014.

%The initialization are the known ones based on NMFInit or the ones of
%Pompili et al, 2014.

% Reference paper: 
% F. Pompili, N. Gillis, P.-A. Absil, F. Glineur,
% "Two algorithms for orthogonal nonnegative matrix 
% factorization with application to clustering".
% Neurocomputing, Vol. 141, pp. 15-25, 2014.

%INPUT:
%pathToSaveResults: Path to save the results of the clustering evaluation.
%numClusters:       Number of Clusters
%tauTV:             Regularization parameter for the TV regularizer, which
%                   is applied as a post processing step.
%alpha0Inp:         Parameter choice for alpha0 (see Pompili, ONP-MF).
%stepBetaInp:       Parameter choice for lsOpts.stepBeta (see Pompili,
%                   ONP-MF).
%iniType:           Initialization Type. String, which is either 'svdbased'
%                   or 'kmeansPP'.
%iniTypePompili:    Specifies how to initialize the matrices within the
%                   code of Pompili et al. Can be a string: Either
%                   'givenUPomp' (given U) or 'givenUVPomp' (given
%                   V)
%indexParam:        Index of the vector of the different sigma values for
%                   the parameter tests
%whichData:         Which kind of data should be used? Should be either
%                   'exampleData' or 'ownData' (see README.md).
%dataPath:          In case of 'ownData', the full path to the *.mat file 
%                   with the dataset and all necessary components is
%                   needed (in case of 'exampleData', just use ''; see
%                   README.md.)

%OUTPUT:
%U:                 Cluster Membership Matrix (of size MxK)
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
    fprintf('Setting the hyperparameters for the NMF...');
    param.iniType                   = iniType; %svdbased, kmeansPP
    param.iniTypePompili            = iniTypePompili;
    param.numClusters               = numClusters;  %Number of Clusters
    param.alpha0Inp                 = alpha0Inp;
    param.stepBetaInp               = stepBetaInp;
%     param.stoppingCrit              = 'relChange'; %Either 'relChange' or 'clusterAssign'
%     param.numEqualInARowMax         = 10; %Maximum number of times that the cluster membership
% %                                  assignments do not change, until the flag is set to
% %                                  true
    param.maxIt                     = 5000;   %Not sure here, Is this the default case for Pompili et al., 2014?
%     param.stopLimUV               = 1e-5;%Stopping Criterion (Relative change) for the matrices U, V.
%     param.rescaling                 = false; % Flag for Rescaling. Default: false
    param.epsStabInit               = [1e-16, 10^35]; %Constant, which is added at the initialization of the NMF
    param.epsLim                    = [1e-16, 10^35]; %Lower and upper projection value during NMF Iterations
    
    %TV Regularization
    param.tauTV = tauTV;
    param.TVFuzzyFlag = true;
%     param.TVDenoiser = 'Driggs';

    fprintf('Done.\n')
    
    %% Preprocessing Step For Datamatrix
    X=projectNegativeAndHighValues(X, param.epsLim);

    %% NMF Initialization
    tStart = tic;
    fprintf(['NMF Initialization: ', param.iniType, '\n']);
    [Uinit,Vinit,warningsUInit,warningsVInit] = NMFInit(X,param.numClusters,size(X,1),size(X,2),...
        param.epsStabInit, param.iniType);
    warnings.UInit = warningsUInit;
    warnings.VInit = warningsVInit;
%     switch param.iniType
%         case 'pompili'
%         case 'pompiliGivenU' %Initialize U with kMeansPP
%             [U,~,warningsUInit,warningsVInit] = NMFInit(X,param.numClusters,size(X,1),size(X,2),...
%                 param.epsStabInit, 'kmeansPP');
%             warnings.UInit = warningsUInit;
%             warnings.VInit = warningsVInit;
%         otherwise
%             [U,V,warningsUInit,warningsVInit] = NMFInit(X,param.numClusters,size(X,1),size(X,2),...
%                 param.epsStabInit, param.iniType);
%             warnings.UInit = warningsUInit;
%             warnings.VInit = warningsVInit;
%     end

    fprintf('Done.\n\n')
    % %% NMF Iterations

%     normFro = @(X1,X2)sqrt(sum(sum((X1-X2).^2)));
%     normFro2 = @(X1,X2)sum(sum((X1-X2).^2)); %#ok<NASGU>

    %% Evaluation Flags

    %Evaluate the cost function at each iteration
%     costFunValueFlag=false;

    %Evaluate Relative Discrepancies
%     relDiscFlag=false;
% 
%     if relDiscFlag
%         relDiscBCTV=zeros(1,param.maxIt);
%     end
% 
%     if costFunValueFlag
%         costFun = @(X,U,V)( 1/2*norm(X-U*V,'fro')^2 ); %#ok<UNRCH>
% 
%         costFunValues=zeros(1,param.maxIt);
%     end

%     %If true, the Feedback of the method shows some hard Clustering results
%     %instead of simply the entries of U.
    hardClustFlag=true; %#ok<NASGU>

    %% ONMF Iterations
%     reverseStr='';
    fprintf('Starting the NMF Iterations based on Pompili et al., 2014 (ONP-MF)\n')
    fprintf('Progress:\n')

    ticTimeAxis=tic;

    switch param.iniTypePompili
        case 'givenUPomp'
            [U,V,costFunValues,numIter]=onpmf_Pompili(X, param.numClusters,...
                param.maxIt, param.alpha0Inp, param.stepBetaInp, Uinit); %#ok<ASGLU>
        case 'givenUVPomp'
            [U,V,costFunValues,numIter]=onpmf_Pompili(X, param.numClusters,...
                param.maxIt, param.alpha0Inp, param.stepBetaInp, Uinit, Vinit); %#ok<ASGLU>
        otherwise
            error('Wrong Choice of param.iniTypePompili!')
    end

    timeAxis=toc(ticTimeAxis);
    fprintf('\n%f minutes needed.\n%d iterations needed.\nFinished!\n\n',...
        round(timeAxis(end)/60,1),numIter)

%     %Compute the (hard) cluster membership matrix U
%     U=indexVectorToClusterMembership(indexVector);

    %% Apply here the TV regularizer and compute clustering quality of U and UTV
    %Be careful here: It has to be a TV postprocessing for hard
    %clusterings!
    fprintf('Applying TV regularization as post processing step...\n')
    UTV = TV_PostProcessing(U, param.TVFuzzyFlag, param.tauTV,...
        data);

    %% Compute Clustering Quality
    clusterQuality = clusteringQuality(U, data);
    clusterQualityTV = clusteringQuality(UTV, data);
%     clusterQuality.relErrorPompili=relErrorPompili;

    msg=['Cluster Quality of U:\n',...
        'Purity: ',num2str(clusterQuality.P), '\nEntropy: ',...
          num2str(clusterQuality.E), '\nVIN: ', num2str(clusterQuality.VIN),...
          '\nVDN: ', num2str(clusterQuality.VDN),'\n'];
        fprintf(msg)

    msg=['Cluster Quality after the TV post processing:',...
        num2str(clusterQualityTV.distUtoIl1), '\nDistance UTV to identity (l2): ',...
          num2str(clusterQualityTV.distUtoIl2), '\nDistance UTV to identity (max. Value): ',...
          num2str(clusterQualityTV.distUtoImax), '\n\nPurity: ',...
          num2str(clusterQualityTV.P), '\nEntropy: ', num2str(clusterQualityTV.E),...
          '\nVIN: ', num2str(clusterQualityTV.VIN), '\nVDN: ',...
          num2str(clusterQualityTV.VDN),'\n'];
    fprintf(msg)

    %% Check Clusterings U and UTV

    warnings.UFuzzy = checkFuzzyClustering(U);
    warnings.UHard = checkHardClustering(convertToHardClustering(U));

    warnings.UTVFuzzy = checkFuzzyClustering(UTV);
    warnings.UTVHard = checkHardClustering(convertToHardClustering(UTV)); %#ok<STRNU>

    %% Compute U^T * U to get the approximative identity matrix
    UU=U'*U;  %#ok<NASGU>
    UUTV = UTV'*UTV; %#ok<NASGU>
%     %% Show Clusters
%     showClusters(U, data, param, true, hardClustFlag, true)

    timeNeeded=toc(tStart);

    %% Save Everything
    data=rmfield(data,'X'); %#ok<NASGU>
    clearvars -except indexParam jj pathToSaveResults timeNeeded numIter fileNameAppend numEqualInARowVec updateChangeVec warnings clusterQualityTV UUTV clusterQualityVec indexParam mfilename param U UTV V W clusterQuality i data stoppingCause timeAxis UU UTV costFunValues
    nameFile=strcat(pathToSaveResults, '\', mfilename, '_test', num2str(indexParam));
    save(nameFile);
    fprintf('Finished.\n')
end