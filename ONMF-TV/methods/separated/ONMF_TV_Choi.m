function [U,UTV,V,clusterQuality,timeNeeded]=ONMF_TV_Choi(pathToSaveResults, numClusters, tauTV, iniType, indexParam, whichData, dataPath)
%Computes the (orthogonal) NMF with orthogonal constaints based on Choi, S.
%Algorithms for orthogonal nonnegative matrix factorization. Proceedings of
%the International Joint Conference on Neural Networks, 2008, pp. 1828–1832.
%The datamatrix X has to be arranged in such a way, that the %observations
%are ordered rowwise.

%INPUT:
%pathToSaveResults: Path to save the results of the clustering evaluation.
%numClusters:       Number of Clusters
%tauTV:             Regularization parameter for the TV regularizer, which
%                   is applied as a post processing step.
%iniType:           Initialization Type. String, which is either 'svdbased'
%                   or 'kmeansPP'.
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
    fprintf('Setting the hyperparameters for the NMF...');
    %Set the main hyper parameters
    param = setMainHyperParam(numClusters);
    
    param.iniType           = iniType; %svdbased, kmeansPP
    param.maxIt             = 600;   %Maximal Iteration Number
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
    [U,V,warningsUInit,warningsVInit] = NMFInit(X,param.numClusters,size(X,1),size(X,2),...
        param.epsStabInit, param.iniType);
    warnings.UInit = warningsUInit;
    warnings.VInit = warningsVInit;

    %Projection Step for multiplicative Update Rules
    U = projectNegativeAndHighValues(U, param.epsStabInit);
    V = projectNegativeAndHighValues(V, param.epsStabInit);

    fprintf('Done.\n\n')
    % %% NMF Iterations

%     normFro = @(X1,X2)sqrt(sum(sum((X1-X2).^2)));
%     normFro2 = @(X1,X2)sum(sum((X1-X2).^2)); %#ok<NASGU>

    %% Evaluation Flags

    %Evaluate the cost function at each iteration
    costFunValueFlag=false;

    %Evaluate Relative Discrepancies
%     relDiscFlag=false;
% 
%     if relDiscFlag
%         relDiscBCTV=zeros(1,param.maxIt);
%     end

    if costFunValueFlag %CAREFUL HERE: Not entirely clear if this is really the cost function of Choi, 2008.
        costFun = @(X,U,V)( 1/2*norm(X-U*V,'fro')^2 ); %#ok<UNRCH>

        costFunValues=zeros(1,param.maxIt);
    end

    %If true, the Feedback of the method shows some hard Clustering results
    %instead of simply the entries of U.
    hardClustFlag=true; %#ok<NASGU>

    %% NMF Iterations
    reverseStr='';
    fprintf('Starting the NMF Iterations based on Choi et al.\n')
    fprintf('Progress:\n')

    i=0;
    stoppingFlag = false;

%     if strcmp(param.stoppingCrit, 'clusterAssign')
    numEqualInARow=0;
%     end

    while ~stoppingFlag

        i=i+1;

%         if relDiscFlag
%             XOld = U*V;
%         end

        U2=U;
        V2=V;

        tic;
        ticTimeAxis=tic;

        %XXXXXXXXXXXXXXXXXXX   Main NMF Iterations %XXXXXXXXXXXXXXXXXXXXXX

        V = V.* ( (U'*X)./( (U'*U)*V ) );

        V=projectNegativeAndHighValues(V, param.epsLim);

        U=U.*( (X*V')./(U*((V*X')*U)) );

        U=projectNegativeAndHighValues(U, param.epsLim);

        %Reduce the projection value with every time step
    %     param.epsLim(1)=1/(20+iii);

        %Project very small or very high values for numerical stability

    %     B = B + param.epsLim(1);

        timeStop=toc(ticTimeAxis);
        if i==1
            timeAxis(i)=timeStop; %#ok<AGROW>
        else
            timeAxis(i)=timeAxis(i-1)+timeStop; %#ok<AGROW>
        end


    %     if param.rescaling
    %         scalingFactors = sqrt(sum(B.*B,1)) + eps;
    %         B = B./repmat(scalingFactors,[size(B,1) 1]);
    %         C = C.*repmat(scalingFactors',[1 size(C,2)]);
    %     end
    %     
        if costFunValueFlag
            costFunValues(i)=costFun(X,U,V); %#ok<UNRCH>
        end
    %     
    %     if relDiscFlag
    %         relDiscBCTV(iii)=normFro(X,XOld)/normFro(XOld, 0); %#ok<*UNRCH>
    %     end


        %Check Stopping Criterion

        [stoppingFlag, updateChange, numEqualInARow, stoppingCause] =...
            stoppingCriterionUV(i,U,U2,V,V2,param,numEqualInARow);

        updateChangeVec(i) = updateChange; %#ok<NASGU,AGROW>
        numEqualInARowVec(i) = numEqualInARow; %#ok<NASGU,AGROW>
%         if strcmp(param.stoppingCrit, 'relChange')
%         [stoppingFlag, updateChange] =...
%             stoppingCriterionUVRelChange(i,U,U2,V,V2,param.stopLimUV,param.maxIt);
%         elseif strcmp(param.stoppingCrit, 'clusterAssign')
%             [stoppingFlag, numEqualInARow]=stoppingCriterionClusterAssign(iii,...
%                 U,U2,numEqualInARow,param.numEqualInARowMax,maxIt);
%         else
%             error('Unknown choice of the stopping criterion.\n')
%         end

        %Feedback of the method
%         reverseStr=methodFeedbackUV(U, i, labels, indexNonzero,...
%             indexGrid, param, updateChange,reverseStr,hardClustFlag,imgMask);
        clusterQuality = clusteringQuality(U, data);
        clusterQualityVec(i) = clusterQuality; %#ok<AGROW>

%         if mod(i, 100)==1
% %             clusterQuality = clusteringQuality(U, labels, indexNonZero);
%             %     close all
% 
%             showClusters(U, data, param, true, hardClustFlag, true)
% 
%     %         figure(1), clf
%     %         mySupTitle=['Columns of Cluster Membership Matrix - Iter: ',...
%     %             num2str(i)];
%     %         suptitle(mySupTitle)
%     %         
%     %         numRowsFig=ceil(param.numClusters/clusterInOneRowFig);
%     % 
%     %         for ii=1:param.numClusters
%     %             subplot(numRowsFig, clusterinOneRowFig, ii)
%     %             imagesc(columnToImage(U(:, ii),indexNonZero, indexGrid));
%     %             axis equal, axis off
%     %             colorbar
%     %             myTitle=['Cluster ', num2str(ii)];
%     %             title(myTitle);
%     %         end
% 
%             %     pause()
%         end
        timeIt=toc;
        timeLeft=round((param.maxIt-i)*timeIt/60,1);
        msg=['Iteration ', num2str(i), '.\nMaximal time left: ', num2str(timeLeft),...
            ' minutes\n\nUpdate Change U: ',...
          num2str(updateChange.U), '\nUpdate Change V: ', num2str(updateChange.V),...
          '\n\nDistance U to identity (l1): ',...
          num2str(clusterQuality.distUtoIl1), '\nDistance U to identity (l2): ',...
          num2str(clusterQuality.distUtoIl2), '\nDistance U to identity (max. Value): ',...
          num2str(clusterQuality.distUtoImax), ...
          '\n\nPurity: ',num2str(clusterQuality.P), '\nEntropy: ',...
          num2str(clusterQuality.E), '\nVIN: ', num2str(clusterQuality.VIN),...
          '\nVDN: ', num2str(clusterQuality.VDN),'\n'];
        fprintf([reverseStr, msg])
        reverseStr=repmat(sprintf('\b'), 1, length(msg)-14);
    end

    fprintf('Main Algorithm Finished!\n')

    %% Perform tasks after the main algorithm
    [UTV,clusterQuality,clusterQualityTV,clusterQualityVec,warnings,UU,...
        UUTV]=tasksAfterMainAlgForSepMethods(stoppingCause,param,U,warnings,...
        data,clusterQualityVec);  %#ok<ASGLU>

    timeNeeded=toc(tStart);

    %% Save Everything
    data=rmfield(data,'X'); %#ok<NASGU>
    clearvars -except indexParam pathToSaveResults timeNeeded fileNameAppend numEqualInARowVec updateChangeVec warnings clusterQualityTV UUTV clusterQualityVec indexParam mfilename param U UTV V W clusterQuality i jj data stoppingCause timeAxis UU costFunValues
    nameFile=strcat(pathToSaveResults, '\', mfilename, '_test', num2str(indexParam));
    save(nameFile);
    fprintf('Finished.\n')
end