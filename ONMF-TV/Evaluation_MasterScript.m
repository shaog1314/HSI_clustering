%% Readme
%This is the masterscript for the evaluation of all considered methods in
%the following, corresponding article:

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

%A short documentation of the code is provided in the file README.md.

%This script is written for general nonnegative hyperspectral datasets. 

%% K-means-TV
fprintf('Evaluation of K-means-TV starts...\n')

%Set the number of Clusters
numClusters = 16; %Default for the example Dataset
%Number of replicates of the internal Matlab K-means algorithm. Default: 1
replicatesKMeans = 1;
%Regularization parameter of the TV postprocessing
tauTVVec = 1;
%Number of Replicates. Default: 30.
numReplicates = 30;
%Which kind of data should be used? Should be either 'exampleData' or 'ownData'.
whichData = 'exampleData';
%In case of 'ownData', the full path to the *.mat file with the dataset and
%all necessary components is needed. (In case of 'exampleData', just use
%'')
dataPath = '';

%Define the location to save the results.
% currentFolder = pwd;
mkdir KMeans_TV
cd KMeans_TV\
pathToSaveResults = pwd;

%Random Seed
rng(42, 'twister');
for i=1:numReplicates
    fprintf('Replicate %d of %d\n', i, numReplicates)
    [U,UTV,V,clusterQuality,timeNeeded]=KMeans_TV(pathToSaveResults, numClusters,...
        replicatesKMeans, tauTVVec, i, whichData, dataPath);
end

%Some Visualization scripts
sepOrCombined = 'separated'; %  Can be either 'separated' for separated methods
%                               or 'combined' for combined methods.
choiceOfMeasure = "Vdn"; %  (see the function printClusterValidationMeasure
%                           for possible cluster validation measures)
printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates,...
    choiceOfMeasure, 1)
plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, 1)

fprintf('Evaluation of K-means-TV ended.\n')

%% ONMF-TV-Choi
fprintf('Evaluation of ONMF_TV_Choi starts...\n')

%Set the number of Clusters
numClusters = 16; %Default for the example Dataset
%Regularization parameter of the TV postprocessing
tauTVVec = 2e-2;
%Number of Replicates. Default: 30.
numReplicates = 30;
%Set initialization procedure
iniType = 'svdbased';
%Which kind of data should be used? Should be either 'exampleData' or 'ownData'.
whichData = 'exampleData';
%In case of 'ownData', the full path to the *.mat file with the dataset and
%all necessary components is needed. (In case of 'exampleData', just use
%'')
dataPath = '';

%Define the location to save the results.
% currentFolder = pwd;
mkdir ONMF_TV_Choi
cd ONMF_TV_Choi\
pathToSaveResults = pwd;

%Random Seed
rng(42, 'twister');
for i=1:numReplicates
    fprintf('Replicate %d of %d\n', i, numReplicates)
    [U,UTV,V,clusterQuality,timeNeeded] = ONMF_TV_Choi(pathToSaveResults,...
        numClusters, tauTVVec, iniType, i, whichData, dataPath);
end

%Some Visualization scripts
sepOrCombined = 'separated'; %  Can be either 'separated' for separated methods
%                               or 'combined' for combined methods.
choiceOfMeasure = "Vdn";%  (see the function printClusterValidationMeasure
%                           for possible cluster validation measures)
printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates,...
    choiceOfMeasure, 1)
plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, 1)

fprintf('Evaluation of ONMF-TV-Choi ended.\n')

%% ONMF-TV-Ding

fprintf('Evaluation of of ONMF-TV-Ding starts...\n')

%Set the number of Clusters
numClusters = 16; %Default for the example Dataset
%Regularization parameter of the TV postprocessing
tauTVVec = 2e-2;
%Number of Replicates. Default: 30.
numReplicates = 30;
%Set initialization procedure
iniType = 'kmeansPP';
%Which kind of data should be used? Should be either 'exampleData' or 'ownData'.
whichData = 'exampleData';
%In case of 'ownData', the full path to the *.mat file with the dataset and
%all necessary components is needed. (In case of 'exampleData', just use
%'')
dataPath = '';


%Define the location to save the results.
% currentFolder = pwd;
mkdir ONMF_TV_Ding
cd ONMF_TV_Ding\
pathToSaveResults = pwd;

%Random Seed
rng(42, 'twister');
for i=1:numReplicates
    fprintf('Replicate %d of %d\n', i, numReplicates)
    [U,UTV,V,clusterQuality,timeNeeded] = ONMF_TV_Ding(pathToSaveResults,...
        numClusters, tauTVVec, iniType, i, whichData, dataPath);
end

%Some Visualization scripts
sepOrCombined = 'separated'; %  Can be either 'separated' for separated methods
%                               or 'combined' for combined methods.
choiceOfMeasure = "Vdn";%  (see the function printClusterValidationMeasure
%                           for possible cluster validation measures)
printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates,...
    choiceOfMeasure, 1)
plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, 1)

fprintf('Evaluation of of ONMF-TV-Ding ended.\n')

%% ONMF-TV-Pompili1

fprintf('Evaluation of ONMF-TV-Pompili1 starts...\n')

%Set the number of Clusters
numClusters = 16; %Default for the example Dataset
%Regularization parameter of the TV postprocessing
tauTVVec = 1;
%Number of Replicates. Default: 30.
numReplicates = 30;
%Set initialization procedure
iniType = 'kmeansPP';
iniTypePompili = 'givenVPomp';
%Which kind of data should be used? Should be either 'exampleData' or 'ownData'.
whichData = 'exampleData';
%In case of 'ownData', the full path to the *.mat file with the dataset and
%all necessary components is needed. (In case of 'exampleData', just use
%'')
dataPath = '';


%Define the location to save the results.
% currentFolder = pwd;
mkdir ONMF_TV_Pompili1
cd ONMF_TV_Pompili1\
pathToSaveResults = pwd;

%Random Seed
rng(42, 'twister');
for i=1:numReplicates
    fprintf('Replicate %d of %d\n', i, numReplicates)
    [U,UTV,V,clusterQuality,timeNeeded] = ONMF_TV_Pompili1(pathToSaveResults,...
        numClusters, tauTVVec, iniType, iniTypePompili, i, whichData, dataPath);
end

%Some Visualization scripts
sepOrCombined = 'separated'; %  Can be either 'separated' for separated methods
%                               or 'combined' for combined methods.
choiceOfMeasure = "Vdn";%  (see the function printClusterValidationMeasure
%                           for possible cluster validation measures)
printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates,...
    choiceOfMeasure, 1)
plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, 1)

fprintf('Evaluation of ONMF-TV-Pompili1 ended.\n')

%% ONMF-TV-Pompili2

fprintf('Evaluation of ONMF-TV-Pompili2 starts...\n')

%Set the number of Clusters
numClusters = 16; %Default for the example Dataset
%Regularization parameter of the TV postprocessing
tauTVVec = 4e-2;
%Number of Replicates. Default: 30.
numReplicates = 30;
%Set initialization procedure
iniType = 'svdbased';
iniTypePompili = 'givenUVPomp';
%Which kind of data should be used? Should be either 'exampleData' or 'ownData'.
whichData = 'exampleData';
%In case of 'ownData', the full path to the *.mat file with the dataset and
%all necessary components is needed. (In case of 'exampleData', just use
%'')
dataPath = '';

%Set further Parameters of the workflow ONMF-TV-Pompili2
alpha0InpVec    = 1e-1;
stepBetaInpVec  = 1.1;

%Define the location to save the results.
% currentFolder = pwd;
mkdir ONMF_TV_Pompili2
cd ONMF_TV_Pompili2\
pathToSaveResults = pwd;


%Random Seed
rng(42, 'twister');
for i=1:numReplicates
    fprintf('Replicate %d of %d\n', i, numReplicates)
    [U,UTV,V,clusterQuality,timeNeeded] = ONMF_TV_Pompili2(pathToSaveResults,...
        numClusters, tauTVVec, alpha0InpVec,...
        stepBetaInpVec, iniType, iniTypePompili, i, whichData, dataPath);
end

%Some Visualization scripts
sepOrCombined = 'separated'; %  Can be either 'separated' for separated methods
%                               or 'combined' for combined methods.
choiceOfMeasure = "Vdn";%  (see the function printClusterValidationMeasure
%                           for possible cluster validation measures)
printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates,...
    choiceOfMeasure, 1)
plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, 1)

fprintf('Evaluation of ONMF-TV-Pompili2 ended.\n')

%% ONMF-TV-Kimura

fprintf('Evaluation of ONMF-TV-Kimura starts...\n')

%Set the number of Clusters
numClusters = 16; %Default for the example Dataset
%Regularization parameter of the TV postprocessing
tauTVVec = 2e-2;
%Number of Replicates. Default: 30.
numReplicates = 30;
%Set initialization procedure
iniType = 'svdbased';
iniTypePompili = 'givenUVPomp';
%Which kind of data should be used? Should be either 'exampleData' or 'ownData'.
whichData = 'exampleData';
%In case of 'ownData', the full path to the *.mat file with the dataset and
%all necessary components is needed. (In case of 'exampleData', just use
%'')
dataPath = '';


%Define the location to save the results.
% currentFolder = pwd;
mkdir ONMF_TV_Kimura
cd ONMF_TV_Kimura\
pathToSaveResults = pwd;

%Random Seed
rng(42, 'twister');
for i=1:numReplicates
    fprintf('Replicate %d of %d\n', i, numReplicates)
    [U,UTV,V,clusterQuality,timeNeeded] = ONMF_TV_Kimura(pathToSaveResults,...
        numClusters, tauTVVec, iniType, i, whichData, dataPath);
end

%Some Visualization scripts
sepOrCombined = 'separated'; %  Can be either 'separated' for separated methods
%                               or 'combined' for combined methods.
choiceOfMeasure = "Vdn";%  (see the function printClusterValidationMeasure
%                           for possible cluster validation measures)
printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates,...
    choiceOfMeasure, 1)
plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, 1)

fprintf('Evaluation of ONMF-TV-Kimura ended.\n')

%% ONMF-TV-Li

fprintf('Evaluation of ONMF-TV-Li starts...\n')

%Set the number of Clusters
numClusters = 16; %Default for the example Dataset
%Regularization parameter of the TV postprocessing
tauTVVec = 3e-2;
%Number of Replicates. Default: 30.
numReplicates = 30;
%Set initialization procedure
iniType = 'svdbased';
%Which kind of data should be used? Should be either 'exampleData' or 'ownData'.
whichData = 'exampleData';
%In case of 'ownData', the full path to the *.mat file with the dataset and
%all necessary components is needed. (In case of 'exampleData', just use
%'')
dataPath = '';

%Define the location to save the results.
% currentFolder = pwd;
mkdir ONMF_TV_Li
cd ONMF_TV_Li\
pathToSaveResults = pwd;

%Set further Parameters of the workflow ONMF-TV-Li
%Computation of the optimal lambda
data=loadData(whichData, dataPath); %load data
X=data.X;
normXFro = norm(X, 'fro');
clear data X
% lambda0 = 0.08;

lambdaVec    = 1e-4 * normXFro;
repeatLines  = 3;

%Random Seed
rng(42, 'twister');
for i=1:numReplicates
    fprintf('Replicate %d of %d\n', i, numReplicates)
    [U,UTV,V,clusterQuality,timeNeeded] = ONMF_TV_Li(pathToSaveResults,...
        numClusters, lambdaVec, repeatLines,...
        tauTVVec, iniType, i, whichData, dataPath);
end

%Some Visualization scripts
sepOrCombined = 'separated'; %  Can be either 'separated' for separated methods
%                               or 'combined' for combined methods.
choiceOfMeasure = "Vdn";%  (see the function printClusterValidationMeasure
%                           for possible cluster validation measures)
printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates,...
    choiceOfMeasure, 1)
plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, 1)

fprintf('Evaluation of ONMF-TV-Li ended.\n')

%% ONMFTV-MUL1

fprintf('Evaluation of ONMFTV-MUL1 starts...\n')

%Set the number of Clusters
numClusters = 16; %Default for the example Dataset
%Regularization parameter of the TV postprocessing
tauTVVec = 5e-3;
%Regularization parameter sigma for the penalty term enforcing the
%orthogonality
sigmaVecFinal=5e-1;
%Number of Replicates. Default: 30.
numReplicates = 30;
%Set initialization procedure
iniType = 'svdbased';
%Which kind of data should be used? Should be either 'exampleData' or 'ownData'.
whichData = 'exampleData';
%In case of 'ownData', the full path to the *.mat file with the dataset and
%all necessary components is needed. (In case of 'exampleData', just use
%'')
dataPath = '';

%Define the location to save the results.
% currentFolder = pwd;
mkdir ONMFTV_MUL1
cd ONMFTV_MUL1\
pathToSaveResults = pwd;

%Random Seed
rng(42, 'twister');
for i=1:numReplicates
    fprintf('Replicate %d of %d\n', i, numReplicates)
    [UTV,V,clusterQuality,timeNeeded] = ONMFTV_MUL1(pathToSaveResults,...
        numClusters, tauTVVec, sigmaVecFinal, sigmaVecFinal, iniType, i,...
        whichData, dataPath);
end

%Some Visualization scripts
sepOrCombined = 'combined'; %  Can be either 'separated' for separated methods
%                               or 'combined' for combined methods.
choiceOfMeasure = "Vdn";%  (see the function printClusterValidationMeasure
%                           for possible cluster validation measures)
printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates,...
    choiceOfMeasure, 1)
plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, 1)

fprintf('Evaluation of ONMFTV-MUL1 ended.\n')

%% ONMFTV-MUL2

fprintf('Evaluation of ONMFTV-MUL2 starts...\n')

%Set the number of Clusters
numClusters = 16; %Default for the example Dataset
%Regularization parameter of the TV postprocessing
tauTVVec = 1e-3;
%Regularization parameter sigma for the penalty term enforcing the
%orthogonality
sigmaVecFinal=1e0;
%Number of Replicates. Default: 30.
numReplicates = 30;
%Set initialization procedure
iniType = 'svdbased';
%Which kind of data should be used? Should be either 'exampleData' or 'ownData'.
whichData = 'exampleData';
%In case of 'ownData', the full path to the *.mat file with the dataset and
%all necessary components is needed. (In case of 'exampleData', just use
%'')
dataPath = '';

%Define the location to save the results.
% currentFolder = pwd;
mkdir ONMFTV_MUL2
cd ONMFTV_MUL2\
pathToSaveResults = pwd;

%Random Seed
rng(42, 'twister');
for i=1:numReplicates
    fprintf('Replicate %d of %d\n', i, numReplicates)
    [UTV,V,clusterQuality,timeNeeded] = ONMFTV_MUL2(pathToSaveResults,...
        numClusters, tauTVVec, sigmaVecFinal, iniType, i,...
        whichData, dataPath);
end

%Some Visualization scripts
sepOrCombined = 'combined'; %  Can be either 'separated' for separated methods
%                               or 'combined' for combined methods.
choiceOfMeasure = "Vdn";%  (see the function printClusterValidationMeasure
%                           for possible cluster validation measures)
printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates,...
    choiceOfMeasure, 1)
plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, 1)

fprintf('Evaluation of ONMFTV-MUL2 ended.\n')

%% ONMFTV-PALM

fprintf('Evaluation of ONMFTV-PALM starts...\n')

%Set the number of Clusters
numClusters = 16; %Default for the example Dataset
%Regularization parameter of the TV postprocessing
tauTVVec = 1e-1;
%Regularization parameter sigma for the penalty term enforcing the
%orthogonality
sigmaVecFinal=1e-1;
%Number of Replicates. Default: 30.
numReplicates = 30;
%Set initialization procedure
iniType = 'svdbased';
%Which kind of data should be used? Should be either 'exampleData' or 'ownData'.
whichData = 'exampleData';
%In case of 'ownData', the full path to the *.mat file with the dataset and
%all necessary components is needed. (In case of 'exampleData', just use
%'')
dataPath = '';


%Define the location to save the results.
% currentFolder = pwd;
mkdir ONMFTV_PALM
cd ONMFTV_PALM\
pathToSaveResults = pwd;

%Random Seed
rng(42, 'twister');
for i=1:numReplicates
    fprintf('Replicate %d of %d\n', i, numReplicates)
    [UTV,V,clusterQuality,timeNeeded] = ONMFTV_PALM(pathToSaveResults, numClusters, tauTVVec,...
        sigmaVecFinal, sigmaVecFinal, iniType, i, whichData, dataPath);
end

%Some Visualization scripts
sepOrCombined = 'combined'; %  Can be either 'separated' for separated methods
%                               or 'combined' for combined methods.
choiceOfMeasure = "Vdn";%  (see the function printClusterValidationMeasure
%                           for possible cluster validation measures)
printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates,...
    choiceOfMeasure, 1)
plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, 1)

fprintf('Evaluation of ONMFTV-PALM ended.\n')

%% ONMFTV-iPALM

fprintf('Evaluation of ONMFTV-iPALM starts...\n')

%Set the number of Clusters
numClusters = 16; %Default for the example Dataset
%Regularization parameter of the TV postprocessing
tauTVVec = 1e-1;
%Regularization parameter sigma for the penalty term enforcing the
%orthogonality
sigmaVecFinal=1e-1;
%Number of Replicates. Default: 30.
numReplicates = 30;
%Set initialization procedure
iniType = 'svdbased';
%Which kind of data should be used? Should be either 'exampleData' or 'ownData'.
whichData = 'exampleData';
%In case of 'ownData', the full path to the *.mat file with the dataset and
%all necessary components is needed. (In case of 'exampleData', just use
%'')
dataPath = '';


%Define the location to save the results.
% currentFolder = pwd;
mkdir ONMFTV_iPALM
cd ONMFTV_iPALM\
pathToSaveResults = pwd;

%Random Seed
rng(42, 'twister');
for i=1:numReplicates
    fprintf('Replicate %d of %d\n', i, numReplicates)
    [UTV,V,clusterQuality,timeNeeded] = ONMFTV_iPALM(pathToSaveResults, numClusters, tauTVVec,...
        sigmaVecFinal, sigmaVecFinal, iniType, i, whichData, dataPath);
end

%Some Visualization scripts
sepOrCombined = 'combined'; %  Can be either 'separated' for separated methods
%                               or 'combined' for combined methods.
choiceOfMeasure = "Vdn";%  (see the function printClusterValidationMeasure
%                           for possible cluster validation measures)
printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates,...
    choiceOfMeasure, 1)
plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, 1)

fprintf('Evaluation of ONMFTV-iPALM ended.\n')

%% ONMFTV-SPRING

fprintf('Evaluation of ONMFTV-SPRING starts...\n')

%Set the number of Clusters
numClusters = 16; %Default for the example Dataset
%Regularization parameter of the TV postprocessing
tauTVVec = 1e-4;
%Regularization parameter sigma for the penalty term enforcing the
%orthogonality
sigmaVecFinal=1e-1;
%Number of Replicates. Default: 30.
numReplicates = 30;
%Set initialization procedure
iniType = 'svdbased';
%Which kind of data should be used? Should be either 'exampleData' or 'ownData'.
whichData = 'exampleData';
%In case of 'ownData', the full path to the *.mat file with the dataset and
%all necessary components is needed. (In case of 'exampleData', just use
%'')
dataPath = '';


%Set further Parameters
subsRatio = 40; %Subsample Ratio. Default Choice: 40
proxTVProj = 1e-3; %Projection of the TV regularization Parameter. Default: 1e-3.

%Define the location to save the results.
% currentFolder = pwd;
mkdir ONMFTV_SPRING
cd ONMFTV_SPRING\
pathToSaveResults = pwd;

%Random Seed
rng(42, 'twister');
for i=1:numReplicates
    fprintf('Replicate %d of %d\n', i, numReplicates)
    [UTV,V,clusterQuality,timeNeeded] = ONMFTV_SPRING(pathToSaveResults, numClusters, tauTVVec,...
        sigmaVecFinal, sigmaVecFinal, subsRatio, proxTVProj, iniType, i,...
        whichData, dataPath);
end

%Some Visualization scripts
sepOrCombined = 'combined'; %  Can be either 'separated' for separated methods
%                               or 'combined' for combined methods.
choiceOfMeasure = "Vdn";%  (see the function printClusterValidationMeasure
%                           for possible cluster validation measures)
printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates,...
    choiceOfMeasure, 1)
plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, 1)

fprintf('Evaluation of ONMFTV-SPRING ended.\n')
