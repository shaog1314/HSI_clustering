function [UTV,V,clusterQuality,timeNeeded]=ONMFTV_iPALM(pathToSaveResults, numClusters, tauTV, sigma1, sigma2, iniType, indexParam, whichData, dataPath)
%Computes the orthogonal NMF with orthogonal and TV penalty terms based 
%on the inertial PALM update rules in Driggs, D.; Tang, J.; Liang, J.;
%Davies, M.; Schönlieb, C.B. SPRING: A fast stochastic proximal
%alternating method for non-smooth non-convex optimization. arXiv preprint
%2020, arXiv: 2002.12266.
%The datamatrix X has to be arranged in such a way, that the
%observations are ordered rowwise.
%This function is primarily designed for parameter tests.

%Reference Papers:
%- Driggs, D.; Tang, J.; Liang, J.;
%  Davies, M.; Schönlieb, C.B. SPRING: A fast stochastic proximal
%  alternating method for non-smooth non-convex optimization. arXiv preprint
%  2020, arXiv: 2002.12266.
%- Pock, T.; Sabach, S. Inertial Proximal Alternating Linearized
%  Minimization (iPALM) for Nonconvex and Nonsmooth Problems. SIAM Journal
%  on Imaging Sciences 2016, 9, 1756–1787.

%The NMF cost function is:
%1/2*||X-U*V||_F^2 + sigma1/2 * ||W^T*U - I||_F^2 + sigma2/2 * ||W-U||_F^2
%+ tauTV/2 * TV(U).

%INPUT:
%pathToSaveResults: Path to save the results of the clustering evaluation.
%numClusters:       Number of Clusters
%tauTV:             Regularization parameter for total Variation
%sigma1:            Regularization parameter for ||I-VW^T||
%sigma2:            Regularization parameter for ||W-V||
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
    param = setMainHyperParam(numClusters);
    param.powerIt           = 5;
%     param.proxTVProj        = proxTVProj;
    param.iniType           = iniType; %svdbased, kmeansPP
    param.maxIt             = 300;
    
    % alpha=5e1;%9e1;
%     param.muB=0;%0
%     param.lambdaB=0;%1e-4;
%     param.muC=1e-1;%1e2
%     param.lambdaC=0;%1e0
    % muX=0;%0
    % lambdaX=0;%0
    param.sigma1=sigma1; %Regularization parameter for ||I-VW^T||
    param.sigma2=sigma2; %Regularization parameter for ||W-V||
    

    %TV Regularization
    param.tauTV=tauTV;%Regularization parameter for total Variation: 8e0
    param.epsilonTV=0;%This has to be zero, since we do not use here the differentiable approximation of the TV penalty
    param.TVWeights=ones(1,param.numClusters);
    
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
    %Initialize W
    W = U; %Consider here maybe some perturbation for W?

    fprintf('Done.\n\n')
    % %% NMF Iterations

%     normFro = @(X1,X2)sqrt(sum(sum((X1-X2).^2)));
%     normFro2 = @(X1,X2)sum(sum((X1-X2).^2)); %#ok<NASGU>

    %% Evaluation Flags

%     %Evaluate the cost function at each iteration
    costFunValueFlag=false;

%     %Evaluate Relative Discrepancies
%     relDiscFlag=false;
% 
%     if relDiscFlag
%         relDiscBCTV=zeros(1,param.maxIt);
%     end

    if costFunValueFlag
        costFun = @(X,U,V,W)(1/2*norm(X-U*V,'fro')^2 + param.sigma1/2 *...
            norm(eye(param.numClusters)-W'*U,'fro')^2 + param.sigma2/2 *...
            norm(W-U,'fro')^2 + param.tauTV* TVPenaltyConv(U, data, param) ); %#ok<UNRCH>
        costFunValues=zeros(1,param.maxIt);
    end

    %If true, the Feedback of the method shows some hard Clustering results
    %instead of simply the entries of U.
    hardClustFlag=true; %#ok<NASGU>


    %% PALM Iterations
    reverseStr='';
    fprintf('Starting the iPALMerror Iterations based on Driggs et al. 2020\n')
    fprintf('Progress:\n')

    i=0;
    stoppingFlag = false;

%     if strcmp(param.stoppingCrit, 'clusterAssign')
        numEqualInARow=0;
%     end

    Uold=U;
    Vold=V;
    Wold=W;


    while ~stoppingFlag%i<202


        %XXXXXXXXXXXXXXXXXXX   Main PALM Iterations %XXXXXXXXXXXXXXXXXXXXXX
        tic;
        ticTimeAxis=tic;


        i=i+1;

        U2=U;
        V2=V;
        W2=W;

        param.iPALM = (i-1)/(i+2);

        %   ----------------Iterations w.r.t. U----------------

%         if relDiscFlag
%             XOld = U*V;
%         end


        Utemp = U + param.iPALM * (U - Uold);
        Uold = Utemp;
        [L_U, ~] = powerMethod(U, V, W, param, 'U');  % estimate Lipschitz constant
        eta_U = 0.9/sum(L_U); %Calculate the step size of the gradient descent step. (The sum: See the derivation...)
        grad_U = Utemp*(V*V') - X*V' + param.sigma1 * (W*(W'*Utemp) - W) +...
            param.sigma2 * (Utemp - W); %Calculate the gradient
        U = Utemp - eta_U*grad_U; %gradient descent step

        %TV Penalty via applying the proximal mapping
        for ii = 1:param.numClusters
            U(:,ii) = proxTV(U(:,ii), data, param.tauTV*eta_U);
        end
        %Nonnegativity Constraint
        U(U<0) = 0;

        %   ----------------Iterations w.r.t. V----------------
        Vtemp = V + param.iPALM * (V - Vold);
        Vold = Vtemp;
        [L_V, ~] = powerMethod(U, V, W, param, 'V');  % estimate Lipschitz constant
        eta_V = 0.9/L_V; %Calculate the step size of the gradient descent step
        grad_V = (U'*U)*Vtemp - U'*X; %Calculate the gradient
        V = Vtemp - eta_V*grad_V; %gradient descent step
        %Nonnegativity Constraint
        V(V<0) = 0;

        %   ----------------Iterations w.r.t. W----------------
        Wtemp = W + param.iPALM * (W - Wold);
        Wold = Wtemp;
        [L_W, ~] = powerMethod(U, V, W, param, 'W');  % estimate Lipschitz constant
        eta_W = 0.9/L_W; %Calculate the step size of the gradient descent step
        grad_W = param.sigma1 *(U*(U'*Wtemp) - U) + param.sigma2*(Wtemp - U); %Calculate the gradient
        W = Wtemp - eta_W*grad_W; %gradient descent step
        %Nonnegativity Constraint
        W(W<0) = 0;

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
            costFunValues(i)=costFun(X,U,V,W); %#ok<UNRCH>
        end
    %     
    %     if relDiscFlag
    %         relDiscBCTV(iii)=normFro(X,XOld)/normFro(XOld, 0); %#ok<*UNRCH>
    %     end


        %Check Stopping Criterion

        [stoppingFlag, updateChange, numEqualInARow, stoppingCause] =...
            stoppingCriterionUVW(i,U,U2,V,V2,W,W2,param,numEqualInARow);

        updateChangeVec(i) = updateChange; %#ok<NASGU,AGROW>
        numEqualInARowVec(i) = numEqualInARow; %#ok<NASGU,AGROW>

%         if strcmp(param.stoppingCrit, 'relChange')
%             [stoppingFlag, updateChange] =...
%                 stoppingCriterionUVWRelChange(i,U,U2,V,V2,W,W2,param.stopLimUVW,param.maxIt);
%         elseif strcmp(param.stoppingCrit, 'clusterAssign')
%             [stoppingFlag, numEqualInARow]=stoppingCriterionClusterAssign(iii,...
%                 U,U2,numEqualInARow,param.numEqualInARowMax,maxIt);
%         else
%             error('Unknown choice of the stopping criterion.\n')
%         end

        %Feedback of the method
%         reverseStr=methodFeedbackUVW(U, i, labels, indexNonzero, indexGrid, param,...
%             updateChange,reverseStr,hardClustFlag,imgMask);

        clusterQuality = clusteringQuality(U, data);
        clusterQualityVec(i) = clusterQuality; %#ok<NASGU,AGROW>

%         if mod(i, 20)==1
% %             clusterQuality = clusteringQuality(U, labels, indexNonzero);
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
          '\nUpdate Change W: ', num2str(updateChange.W), '\n\nDistance U to identity (l1): ',...
          num2str(clusterQuality.distUtoIl1), '\nDistance U to identity (l2): ',...
          num2str(clusterQuality.distUtoIl2), '\nDistance U to identity (max. Value): ',...
          num2str(clusterQuality.distUtoImax), '\n\nPurity: ',...
          num2str(clusterQuality.P), '\nEntropy: ', num2str(clusterQuality.E),...
          '\nVIN: ', num2str(clusterQuality.VIN), '\nVDN: ',...
          num2str(clusterQuality.VDN),'\n'];
        fprintf([reverseStr, msg])
        reverseStr=repmat(sprintf('\b'), 1, length(msg)-15);
    end

    fprintf('\nFinished!\n')
    switch stoppingCause
        case 'relChange'
            fprintf('Stopping Cause: Relative Change of the matrices U, V and W.\n')
        case 'clusterAssign'
            fprintf('Stopping Cause: No change of cluster assignments during %d iterations.\n',...
                param.numEqualInARowMax)
        case 'maxIt'
            fprintf('Stopping Cause: Maximum number of iterations (%d) reached.\n',param.maxIt)
    end

    UTV = U;

    %% Compute Clustering Quality
    clusterQualityTV = clusteringQuality(UTV, data); %#ok<NASGU>

        %% Check Clusterings U and UTV

    warnings.UTVFuzzy = checkFuzzyClustering(UTV);
    warnings.UTVHard = checkHardClustering(convertToHardClustering(UTV)); %#ok<STRNU>

    %% Compute final matrices
    UU=U'*U;  %#ok<NASGU>

    UUTV = UTV'*UTV; %#ok<NASGU>

    timeNeeded=toc(tStart);

    %% Save Everything
    data=rmfield(data,'X'); %#ok<NASGU>
    clearvars -except indexParam jjj pathToSaveResults timeNeeded fileNameAppend numEqualInARowVec updateChangeVec warnings clusterQualityTV UUTV clusterQualityVec indexParam mfilename param U UTV V W clusterQuality i data stoppingCause timeAxis UU costFunValues
    nameFile=strcat(pathToSaveResults, '\', mfilename, '_test', num2str(indexParam));
    save(nameFile);
    fprintf('Finished.\n')
end

function [columnUOut]   =   proxTV(columnUIn, data, tauTV)

    param.verbose = 0;
    imageU = columnToImage(columnUIn, data);
    [imageUTV,~] = prox_tv_19(imageU, tauTV, param);
    columnUOut = imageToColumn(imageUTV, data);
    
end
