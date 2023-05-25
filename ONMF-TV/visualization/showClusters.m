function showClusters(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, whichReplicate)
%Extracts and prints the clusterings of the best performed replicates with
%repsect to the specified cluster validation measure from an evaluation of
%one of the ONMF models.

%INPUT:
%pathToSaveResults: Full path of the results
%sepOrCombined:     Can be either 'separated' for separated methods or
%                   'combined' for combined methods.
%numReplicates:     Number of replicates
%choiceOfMeasure:   String. Specifies the measure, which shall
%                   be plotted.
%                   Expected strings:  "E","P","VI","VD","Vin","Vdn",
%                   "distUtoIl1","distUtoIl2","distUtoImax".
%                   E: Entropy
%                   P: Purity
%                   VI: Variation of Information
%                   VD: Van Dongen criterion
%                   Vin: Normalized Variation of Information
%                   Vdn: Normalized van Dongen criterion
%whichReplicate:    Integer. Choice of replicate, of which the cluster
%                   validation measure should be printed (optional)

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.
    
    clusterQualityMat = getClusterValidationMeasure(pathToSaveResults,...
            sepOrCombined, numReplicates, choiceOfMeasure);
        
    baseNameOfFiles = getBaseNameOfFiles(pathToSaveResults);
    
    if strcmp(choiceOfMeasure, 'E')
        selectedMeasure = 'Entropy (E)';
    elseif strcmp(choiceOfMeasure, 'VI')
        selectedMeasure = 'Variation of Information (VI)';
    elseif strcmp(choiceOfMeasure, 'VD')
        selectedMeasure = 'Van Dongen criterion (VD)';
    elseif strcmp(choiceOfMeasure, 'Vin')
        selectedMeasure = 'Normalized Variation of Information (Vin)';
    elseif strcmp(choiceOfMeasure, 'Vdn')
        selectedMeasure = 'Normalized van Dongen criterion (Vdn)';
    else
        error('Cluster Validation measure not recognized')
    end
    
    figure
    
    if strcmp(sepOrCombined, 'combined')
        numSubFigs = 2;
    elseif strcmp(sepOrCombined, 'separated')
        numSubFigs = 3;
    end
    
    if nargin == 5
        numSubFigs = numSubFigs + 1;
    end
    
    %Plot ground truth
    subPlotGT = subplot(1, numSubFigs, 1);
    fileName=strcat(baseNameOfFiles, num2str(1), '.mat');
    load(fileName, 'data')
    subPlotGroundTruth = imagesc(data.labels);

    numClasses = max(max(data.labels));
    classlabelNames = cell(1,numClasses);
    for k = 1:numClasses
        labelName = ['Class ', num2str(k)];
        classlabelNames{k} = labelName;
    end

    set(subPlotGroundTruth, 'AlphaData', data.imgMask.labels);
    axis equal, axis off

    %Pink: Dark to Bright
%     mapOfColor = [58, 19, 83;...
%         88, 29, 127;...
%         135, 46, 147;...
%         200, 72, 138;...
%         235, 117, 144;...
%         246, 181, 164]./255;
    myMap = parula(numClasses);
    colormap(subPlotGT, myMap)

    %Adjust colorbar
    startingPoint = 1 + (numClasses - 1) / (2*numClasses);
    slidingInterval = (numClasses - 1) / (numClasses);

    ticksVec = zeros(1, numClasses);
    ticksVec(1) = startingPoint;
    for k=2:numClasses
        ticksVec(k) = startingPoint + slidingInterval*(k-1);
    end
    caxis([1 numClasses])

    colorbar('Ticks', ticksVec, 'TickLabels', classlabelNames,...
        'Location', 'westoutside', 'FontSize', 14);
%     colorbar('Location', 'westoutside', 'FontSize', 14);

    title('Ground Truth', 'FontSize', 16)
    
    
    if strcmp(sepOrCombined, 'combined')
        %Plot Clustering of best Replicate w.r.t. specified cluster
        %validation measure
        subPlotComb = subplot(1, numSubFigs, 2);
        [~, minInd] = min(clusterQualityMat);
        fileName=strcat(baseNameOfFiles, num2str(minInd), '.mat');
        load(fileName, 'UTV', 'data')
        
        %Sort the columns of U
        normUVec=sum(UTV.^2, 1);
        [~, sortingIndex] = sort(normUVec, 'descend');
        Usorted = UTV(:,sortingIndex);
        
        %Convert to an image
        indexVector=clusterMembershipToIndexVector(Usorted);
        hardClusteringImage = columnToImage(indexVector, data);
        subPlotClusteringTV=imagesc(hardClusteringImage);
        set(subPlotClusteringTV, 'AlphaData', data.imgMask.labels);
        
        axis equal, axis off
        
        %Determine number of Clusters
        maxNumClusters = max(max(hardClusteringImage));
        clusterlabelNames = cell(1,maxNumClusters);
        for k = 1:maxNumClusters
            labelName = ['Cluster ', num2str(k)];
            clusterlabelNames{k} = labelName;
        end
        
        myMap = parula(maxNumClusters);
        colormap(subPlotComb, myMap)
        
        %Adjust colorbar
        startingPoint = 1 + (maxNumClusters - 1) / (2*maxNumClusters);
        slidingInterval = (maxNumClusters - 1) / (maxNumClusters);

        ticksVec = zeros(1, maxNumClusters);
        ticksVec(1) = startingPoint;
        for k=2:maxNumClusters
            ticksVec(k) = startingPoint + slidingInterval*(k-1);
        end
        caxis([1 maxNumClusters])
        
        colorbar('Ticks', ticksVec, 'TickLabels', clusterlabelNames,...
            'Location', 'westoutside', 'FontSize', 14);
        title(['Best Replicate (Test ', int2str(minInd), ')'], 'FontSize', 16)
        
    elseif strcmp(sepOrCombined, 'separated')
        %Plot Clustering of best Replicate w.r.t. specified cluster
        %validation measure
        subPlotSep1 = subplot(1, numSubFigs, 2);
        [~, minInd] = min(clusterQualityMat(:,1));
        fileName=strcat(baseNameOfFiles, num2str(minInd), '.mat');
        load(fileName, 'U', 'data')
        
        %Sort the columns of U
        normUVec=sum(U.^2, 1);
        [~, sortingIndex] = sort(normUVec, 'descend');
        Usorted = U(:,sortingIndex);
        
        %Convert to an image
        indexVector=clusterMembershipToIndexVector(Usorted);
        hardClusteringImage = columnToImage(indexVector, data);
        subPlotClustering=imagesc(hardClusteringImage);
        set(subPlotClustering, 'AlphaData', data.imgMask.labels);
        
        axis equal, axis off
        
        %Determine number of Clusters
        maxNumClusters = max(max(hardClusteringImage));
        clusterlabelNames = cell(1,maxNumClusters);
        for k = 1:maxNumClusters
            labelName = ['Cluster ', num2str(k)];
            clusterlabelNames{k} = labelName;
        end
        
        myMap = parula(maxNumClusters);
        colormap(subPlotSep1, myMap)
        
        %Adjust colorbar
        startingPoint = 1 + (maxNumClusters - 1) / (2*maxNumClusters);
        slidingInterval = (maxNumClusters - 1) / (maxNumClusters);

        ticksVec = zeros(1, maxNumClusters);
        ticksVec(1) = startingPoint;
        for k=2:maxNumClusters
            ticksVec(k) = startingPoint + slidingInterval*(k-1);
        end
        caxis([1 maxNumClusters])
        
        colorbar('Ticks', ticksVec, 'TickLabels', clusterlabelNames,...
            'Location', 'westoutside', 'FontSize', 14);
        title(['Best Replicate without TV Postproc. (Test ', int2str(minInd), ')'], 'FontSize', 16)
        
        %Plot Clustering of best Replicate w.r.t. specified cluster
        %validation measure
        subPlotSep2 = subplot(1, numSubFigs, 3);
        [~, minInd] = min(clusterQualityMat(:,2));
        fileName=strcat(baseNameOfFiles, num2str(minInd), '.mat');
        load(fileName, 'UTV', 'data')
        
        %Sort the columns of U
        normUVec=sum(UTV.^2, 1);
        [~, sortingIndex] = sort(normUVec, 'descend');
        Usorted = UTV(:,sortingIndex);
        
        %Convert to an image
        indexVector=clusterMembershipToIndexVector(Usorted);
        hardClusteringImage = columnToImage(indexVector, data);
        subPlotClusteringTV=imagesc(hardClusteringImage);
        set(subPlotClusteringTV, 'AlphaData', data.imgMask.labels);
        
        axis equal, axis off
        
        %Determine number of Clusters
        maxNumClusters = max(max(hardClusteringImage));
        clusterlabelNames = cell(1,maxNumClusters);
        for k = 1:maxNumClusters
            labelName = ['Cluster ', num2str(k)];
            clusterlabelNames{k} = labelName;
        end
        
        myMap = parula(maxNumClusters);
        colormap(subPlotSep2, myMap)
        
        %Adjust colorbar
        startingPoint = 1 + (maxNumClusters - 1) / (2*maxNumClusters);
        slidingInterval = (maxNumClusters - 1) / (maxNumClusters);

        ticksVec = zeros(1, maxNumClusters);
        ticksVec(1) = startingPoint;
        for k=2:maxNumClusters
            ticksVec(k) = startingPoint + slidingInterval*(k-1);
        end
        caxis([1 maxNumClusters])
        
        colorbar('Ticks', ticksVec, 'TickLabels', clusterlabelNames,...
            'Location', 'westoutside', 'FontSize', 14);
        title(['Best Replicate with TV Postproc. (Test ', int2str(minInd), ')'], 'FontSize', 16)
    end
    
    %Plot Clustering of specified Replicate (optional)
    if nargin == 5
        if strcmp(sepOrCombined, 'combined')
            subPlotSepSpec = subplot(1, numSubFigs, 3);
        elseif strcmp(sepOrCombined, 'separated')
            subPlotSepSpec = subplot(1, numSubFigs, 4);
        end
        fileName=strcat(baseNameOfFiles, num2str(whichReplicate), '.mat');
        load(fileName, 'U', 'data')

        %Sort the columns of U
        normUVec=sum(U.^2, 1);
        [~, sortingIndex] = sort(normUVec, 'descend');
        Usorted = U(:,sortingIndex);

        %Convert to an image
        indexVector=clusterMembershipToIndexVector(Usorted);
        hardClusteringImage = columnToImage(indexVector, data);
        subPlotClusteringSpec=imagesc(hardClusteringImage);
        set(subPlotClusteringSpec, 'AlphaData', data.imgMask.labels);

        axis equal, axis off

        %Determine number of Clusters
        maxNumClusters = max(max(hardClusteringImage));
        clusterlabelNames = cell(1,maxNumClusters);
        for k = 1:maxNumClusters
            labelName = ['Cluster ', num2str(k)];
            clusterlabelNames{k} = labelName;
        end
        
        myMap = parula(maxNumClusters);
        colormap(subPlotSepSpec, myMap)

        %Adjust colorbar
        startingPoint = 1 + (maxNumClusters - 1) / (2*maxNumClusters);
        slidingInterval = (maxNumClusters - 1) / (maxNumClusters);

        ticksVec = zeros(1, maxNumClusters);
        ticksVec(1) = startingPoint;
        for k=2:maxNumClusters
            ticksVec(k) = startingPoint + slidingInterval*(k-1);
        end
        caxis([1 maxNumClusters])

        colorbar('Ticks', ticksVec, 'TickLabels', clusterlabelNames,...
            'Location', 'westoutside', 'FontSize', 14);
        title(['Test ', int2str(whichReplicate)], 'FontSize', 16)
    end
%     
%     %Add suptitle
%     if strcmp(sepOrCombined, 'combined')
%         suptitle('Comparison of clusterings of a combined ONMF model')
%     elseif strcmp(sepOrCombined, 'separated')
%         suptitle('Comparison of clusterings of a separated ONMF model')
%     end
end
