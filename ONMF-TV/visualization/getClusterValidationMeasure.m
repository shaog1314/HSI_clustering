function clusterQualityMat = getClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
%Extracts the specified cluster validation measure from an evaluation of
%one of the ONMF models.

%INPUT:
%pathToSaveResults: Full path of the results
%sepOrCombined:     Can be either 'separated' for separated methods or
%                   'combined' for combined methods.
%numReplicates:     Number of replicates
%choiceOfMeasure:   String. Specifies the measure, which shall
%                   be plotted.
%                   Expected strings:  "E","VI","VD","Vin","Vdn"
%                   E: Entropy
%                   VI: Variation of Information
%                   VD: van Dongen criterion
%                   Vin: Normalized Variation of Information
%                   Vdn: Normalized van Dongen criterion

%OUTPUT:
%clusterQualityMat: Vector of the cluster validation measures (combined
%                   methods) or matrix with the first column being the 
%                   results without TV postprocessing and the second column 
%                   being the results with TV postprocessing.

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.
    
    %Define the available quality measures as strings and choose the
    %specified ones (see choiceOfMeasure in the input)
    availableMeasures = ["E","P","VI","VD","Vin","Vdn","distUtoIl1",...
        "distUtoIl2","distUtoImax"];
    measureIndex = find(contains(availableMeasures,choiceOfMeasure));
%     numMeasureInd = length(measureIndex);
    
%     pathToExp = getPathToExp(datasetName, methodName, whichTest);
%     cd(pathToExp)
    
%     infoOfElements = dir('*test*.mat');
%     numExp=numel(infoOfElements);
%     baseNameOfFiles = infoOfElements(1).name(1:end-5); %Take always the first 
    %file and erase the last 5 characters to get the base name
    
%     fileName=strcat(baseNameOfFiles, num2str(1), '.mat');
%     load(fileName,'clusterQualityTV');
    if strcmp(sepOrCombined, 'combined')
        clusterQualityMat = zeros(numReplicates, 1);
    elseif strcmp(sepOrCombined, 'separated')
        clusterQualityMat = zeros(numReplicates, 2);
    end
    
    if strcmp(sepOrCombined, 'combined')
        %Load the cluster Quality from all experiments
%         clusterQualityVecsInCell = cell(1,length(numExp));
        for ii=1:numReplicates
            baseNameOfFiles = getBaseNameOfFiles(pathToSaveResults);
            fileName=strcat(baseNameOfFiles, num2str(ii), '.mat');
            load(fileName,'clusterQualityTV');
%             clusterQualityInCell = 
            clusterQualityCellTV = struct2cell(clusterQualityTV);
            clusterQualityMat(ii,1) = clusterQualityCellTV{measureIndex};
%             clusterQualityVecsInCell{ii} = updateChangeVec;
%             clusterQualityCell=struct2cell(clusterQualityVecsInCell{ii});
        end
    elseif strcmp(sepOrCombined, 'separated')
        %Load the cluster Quality from all experiments
%         clusterQualityVecsInCell = cell(1,length(numExp));
        for ii=1:numReplicates
            baseNameOfFiles = getBaseNameOfFiles(pathToSaveResults);
            fileName=strcat(baseNameOfFiles, num2str(ii), '.mat');
            load(fileName,'clusterQualityTV');
            load(fileName,'clusterQuality');
%             clusterQualityInCell = 
            clusterQualityCellTV = struct2cell(clusterQualityTV);
            clusterQualityCell = struct2cell(clusterQuality);
            clusterQualityMat(ii,1) = clusterQualityCell{measureIndex};
            clusterQualityMat(ii,2) = clusterQualityCellTV{measureIndex};
%             clusterQualityVecsInCell{ii} = updateChangeVec;
%             clusterQualityCell=struct2cell(clusterQualityVecsInCell{ii});
        end
    end
%     meanQualityVec = mean(clusterQualityMat);
    
end