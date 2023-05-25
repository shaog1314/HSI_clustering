function printClusterValidationMeasure(pathToSaveResults, sepOrCombined, numReplicates, choiceOfMeasure, whichReplicate)
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
%                   VD: Van Dongen criterion
%                   Vin: Normalized Variation of Information
%                   Vdn: Normalized van Dongen criterion
%whichReplicate:    Integer. Choice of replicate, of which the cluster
%                   validation measure should be printed (optional)


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

    fprintf('<strong>Selected Measure: %s</strong>\n', selectedMeasure);

    
    
    if strcmp(sepOrCombined, 'combined')
        [~, minInd] = min(clusterQualityMat);
        if nargin == 4
            msg=['Best Replicate (Test ', num2str(minInd) '): ', num2str(min(clusterQualityMat)), '\n', 'Mean Value: ',...
                num2str(mean(clusterQualityMat)), '\n\n'];
        elseif nargin > 4
            msg=['Best Replicate (Test ', num2str(minInd) '): ', num2str(min(clusterQualityMat)), '\n', 'Mean Value: ',...
                num2str(mean(clusterQualityMat)), '\n', 'Test ', num2str(whichReplicate),...
                ': ', num2str(clusterQualityMat(whichReplicate)), '\n\n'];
        end
        fprintf(msg)
    elseif strcmp(sepOrCombined, 'separated')
        [~, minInd] = min(clusterQualityMat);
        if nargin == 4
            msg=['Results without TV Postprocessing:\n', 'Best Replicate (Test ', num2str(minInd(1)) '): ',...
                num2str(min(clusterQualityMat(:,1))), '\n', 'Mean Value: ',...
                num2str(mean(clusterQualityMat(:,1))), '\n\n',...
                'Results with TV Postprocessing:\n', 'Best Replicate (Test ', num2str(minInd(2)) '): ',...
                num2str(min(clusterQualityMat(:,2))), '\n', 'Mean Value: ',...
                num2str(mean(clusterQualityMat(:,2))), '\n\n'];
        elseif nargin > 4
            msg=['Results without TV Postprocessing:\n', 'Best Replicate (Test ', num2str(minInd(1)) '): ',...
                num2str(min(clusterQualityMat(:,1))), '\n', 'Mean Value: ',...
                num2str(mean(clusterQualityMat(:,1))), '\n', 'Test ', num2str(whichReplicate),...
                ': ', num2str(clusterQualityMat(whichReplicate,1)), '\n\n',...
                'Results with TV Postprocessing:\n', 'Best Replicate (Test ', num2str(minInd(2)) '): ',...
                num2str(min(clusterQualityMat(:,2))), '\n', 'Mean Value: ',...
                num2str(mean(clusterQualityMat(:,2))), '\n', 'Test ', num2str(whichReplicate),...
                ': ', num2str(clusterQualityMat(whichReplicate,2)), '\n\n'];
        end
        fprintf(msg)
    end
    
end