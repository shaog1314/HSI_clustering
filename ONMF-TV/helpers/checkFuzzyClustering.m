function warningsFound=checkFuzzyClustering(U)
%Check if all the datapoints in U are assigned to unique centroids and give
%a warning to the command window if this is not the case.

%INPUT:
%U:                 Cluster Membership Matrix

%OUTPUT:
%warningsFound:     True, if warnings were found, false otherwise.

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.


    numEqualAssign=sum(U==max(U,[],2),2);
%     warningsFound = false;
    
    warningsFound.neg = false;
    warningsFound.noDat = false;
    warningsFound.prob = false;
    
    if ~isempty(U(U<0))
        warning('Negative entries found in U!\n')
        warningsFound.neg = true;
    end
    if any(all(~U,1))
        warning('Some cluster do not contain any datapoint!\n')
        warningsFound.noDat = true;
    end
    if ~isempty(numEqualAssign(numEqualAssign>1))
        warning('Some datapoints in U have the same probability assigned to different centroids.\n')
        warningsFound.prob = true;
    end
%     else
%         cprintf('green', 'No warnings found for U!\n')
%     end
    
end