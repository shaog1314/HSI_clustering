function [UNew, VNew, warningsUInit, warningsVInit]=checkClusteringInit(U, V)
%Checks several properties of the obtained cluster membership matrix and
%the matrix V, which contains row-wise the centroids. It checks the fuzzy
%and the hard clustering properties and do some suitable projection steps
%if needed.

%INPUT:
%U:                 Cluster Membership Matrix
%V:                 Matrix, which contains the centroids row-wise.

%OUTPUT:
%UNew:              New Cluster Membership Matrix
%VNew:              New Matrix, which contains the centroids row-wise.
%warningsUInit:     True, if warnings have been found for U
%warningsVInit:     True, if warnings have been found for V

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    %Check Fuzzy Clustering
    
    numEqualAssign=sum(U==max(U,[],2),2);
%     warningsUInit = false;
    
    warningsUInit.neg = false;
    warningsUInit.noDat = false;
    warningsUInit.prob = false;
    warningsUInit.multNonZEntries = false;
%     warningsVInit = false;
    
    warningsVInit.neg = false;
    
    UNew = U;
    VNew = V;
    
    if ~isempty(U(U<0))
        warning('At Initialization: Negative entries found in U! They will be set to zero.\n')
        UNew(UNew<0) = 0;
        warningsUInit.neg = true;
    end
    if any(all(~U,1))
        warning('At Initialization: Some cluster do not contain any datapoint!\n')
        warningsUInit.noDat = true;
    end
    if ~isempty(numEqualAssign(numEqualAssign>1))
        warning('At Initialization: Some datapoints in U have the same probability aignment to different centroids.\n')
        warningsUInit.prob = true;
    end
%     else
%         cprintf('green', 'No warnings found for the initialization of the fuzzy Clustering!\n')
%     end
    
    %Check Hard Clustering
    Uhard=convertToHardClustering(U);
    
    numEqualAssign=sum(Uhard==max(Uhard,[],2),2);

    if ~isempty(Uhard(Uhard<0))
        warning('At Initialization: Negative entries found in the hard clustering!\n')
        warningsUInit.neg = true;
    end
    if any(all(~Uhard,1))
        warning('At Initialization: Some cluster do not contain any datapoint in the hard clustering!\n')
        warningsUInit.noDat = true;
    end
    if ~isempty(numEqualAssign(numEqualAssign>1))
        warning('At Initialization: Some datapoints in the hard clustering have the same probability aignment to different centroids.\n')
        warningsUInit.prob = true;
    end
    if any(sum(~~Uhard,2)>1)
        warning('There are multiple nonzero entries in the rows of the hard clustering matrix U.')
        warningsUInit.multNonZEntries = true;
    end
%     else
%         cprintf('green', 'No warnings found for the hard clustering at initialization!\n')
%     end
    
    %Check V
    if any(any(V<0))
        warning('The Matrix V contains at initialization negative elements! They will be set to zero.\n')
        VNew(VNew<0) = 0;
        warningsVInit.neg = true;
    else
        fprintf('No warnings found for the matrix V at initialization!\n')
    end
    
%     UNew=U;
%     VNew=V;
    %Check V
%     
%     if any(any(V==0))
%         warning('The Matrix V contains at initialization zero elements. A suitable projection step will be performed to avoid zero entries.\n')
%         VNew=projectNegativeAndHighValues(V, epsStabInit);
%     else
%         cprintf('green', 'No warnings found for the matrix V at initialization!\n')
%         VNew=V;
%     end
    
end