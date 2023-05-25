function baseNameOfFiles = getBaseNameOfFiles(pathToExp)
%Computes the path to the experiments for a given name of the dataset and
%the name of the method.

%INPUT:
%datasetName:       Name of dataset to define the path of the experiments
%methodName:        Choice of method.

%OUTPUT:
%pathToExp:         Path to the Experiment

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.
    
%Check Input
    if nargin ~= 1
        error('Wrong number of input arguments')
    end
    
    cd(pathToExp)
    infoOfElements = dir('*test*.mat');
%     numExp=numel(infoOfElements);
    baseNameOfFiles = infoOfElements(1).name; 
    baseNameOfFiles = baseNameOfFiles(1:end-5);%Take always the first 
    %file and erase the last 5 characters to get the base name
end