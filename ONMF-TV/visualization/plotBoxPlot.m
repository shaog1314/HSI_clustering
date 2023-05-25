function plotBoxPlot(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure)
%Plots the box plot of all replicates from an evaluation of
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
%                   VD: Van Dongen criterion
%                   Vin: Normalized Variation of Information
%                   Vdn: Normalized van Dongen criterion

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    
    clusterQualityMat = getClusterValidationMeasure(pathToSaveResults,...
        sepOrCombined, numReplicates, choiceOfMeasure);
    
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
    if strcmp(sepOrCombined, 'separated')
        boxplot(clusterQualityMat, 'BoxStyle', 'outline', 'Labels', {'Without TV','With TV'})
    elseif strcmp(sepOrCombined, 'combined')
        boxplot(clusterQualityMat, 'BoxStyle', 'outline', 'Labels', {'Combined Method (with TV)'})
    end
    ylabel(selectedMeasure)
    grid on
%     
%     %Get the performance matrices. Matrices are  correct; I've checked
%     %them.
%     [performanceMatrixSep, performanceMatrixComb] =...
%         getPerfMatOfFinalEval(datasetName, numExp, whichTest,...
%         choiceOfMeasures);
%     
%     
%     
%     figure
%     boxplot(performanceMatrixComb, 'BoxStyle', 'outline')
%     grid on
%     figure
%     boxplot(performanceMatrixSep(:,:,1), 'BoxStyle', 'outline')
%     grid on
%     figure
%     boxplot(performanceMatrixSep(:,:,2), 'BoxStyle', 'outline')
%     grid on
%     hold on
    
%     figure
%     notBoxPlot(performanceMatrixComb,'style','line');
%     grid on
    
%               Examples
        
% %             Example 1: Basic grouped box plot with legend
%         
%               y = randn(50,3,3);
%               x = [1 2 3.5];
%               y(1:25) = NaN;
%         
%               figure;
%               h = iosr.statistics.boxPlot(x,y,...
%                   'symbolColor','k',...
%                   'medianColor','k',...
%                   'symbolMarker',{'+','o','d'},...
%                   'boxcolor',{[1 0 0]; [0 1 0]; [0 0 1]},...
%                   'groupLabels',{'y1','y2','y3'},...
%                   'showLegend',true);
%               box on
%         
%             Example 2: Grouped box plot with overlayed data
        
%               figure;
%               iosr.statistics.boxPlot(x,y,...
%                   'symbolColor','k',...
%                   'medianColor','k',...
%                   'symbolMarker',{'+','o','d'},...
%                   'boxcolor','auto',...
%                   'showScatter',true);
%               box on
%         
% %             Example 3: Grouped box plot with displayed sample sizes
% %               and variable widths
%         
%               figure;
%               iosr.statistics.boxPlot(x,y,...
%                   'medianColor','k',...
%                   'symbolMarker',{'+','o','d'},...
%                   'boxcolor','auto',...
%                   'sampleSize',true,...
%                   'scaleWidth',true);
%               box on
%         
% %             Example 4: Grouped notched box plot with x separators and
% %               hierarchical labels
%         
%               figure;
%               iosr.statistics.boxPlot({'A','B','C'},y,...
%                   'notch',true,...
%                   'medianColor','k',...
%                   'symbolMarker',{'+','o','d'},...
%                   'boxcolor','auto',...
%                   'style','hierarchy',...
%                   'xSeparator',true,...
%                   'groupLabels',{{'Group 1','Group 2','Group 3'}});
%               box on
%         
% %             Example 5: Box plot with legend labels from data
%         
%               % load data
%               % (requires Statistics or Machine Learning Toolbox)
%               load carbig
%         
%               % arrange data
%               [y,x,g] = iosr.statistics.tab2box(Cylinders,MPG,when);
%           
%               % sort
%               IX = [1 3 2]; % order
%               g = g{1}(IX);
%               y = y(:,:,IX);
%         
%               % plot
%               figure
%               h = iosr.statistics.boxPlot(x,y,...
%                   'symbolColor','k','medianColor','k','symbolMarker','+',...
%                   'boxcolor',{[1 1 1],[.75 .75 .75],[.5 .5 .5]},...
%                   'scalewidth',true,'xseparator',true,...
%                   'groupLabels',g,'showLegend',true);
%               box on
%               title('MPG by number of cylinders and period')
%               xlabel('Number of cylinders')
%               ylabel('MPG')
%         
% %             Example 6: Box plot calculated from weighted quantiles
%         
%               % load data
%               load carbig
%               
%               % random weights
%               weights = rand(size(MPG));
%               
%               % arrange data
%               [y,x,g] = iosr.statistics.tab2box(Cylinders,MPG,when);
%               weights_boxed = iosr.statistics.tab2box(Cylinders,weights,when);
%               
%               % plot
%               figure
%               h = iosr.statistics.boxPlot(x,y,'weights',weights_boxed);
%         
% %             Example 7: Draw a violin plot
%               y = randn(50,3,3);
%               x = [1 2 3.5];
%               y(1:25) = NaN;
%               figure('color','w');
%               h2 = iosr.statistics.boxPlot(x,y, 'showViolin', true, 'boxWidth', 0.025, 'showOutliers', false);
%               box on
%     
%     figure
%     plot(1:numExp,performanceMatrixComb,'LineWidth', 2);
%     legend
%     grid on
% 
%     j=7;
%     
%     figure
%     plot(1:numExp,performanceMatrixSep(:,j,1), 'LineWidth', 2);
%     hold on
%     plot(1:numExp,performanceMatrixSep(:,j,2), 'LineWidth', 2);
%     legend
%     grid on


end
