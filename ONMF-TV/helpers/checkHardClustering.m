function warningsFound=checkHardClustering(Uhard)
%Check if all the datapoints of the hardclustering obtained based on U are
%assigned to unique centroids and give a warning to the command window if
%this is not the case.

%INPUT:
%Uhard:             Cluster Membership Matrix

%OUTPUT:
%warningsFound:     True, if warnings were found, false otherwise.

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.


    numEqualAssign=sum(Uhard==max(Uhard,[],2),2);
%     warningsFound = false;
    
    warningsFound.neg = false;
    warningsFound.noDat = false;
    warningsFound.prob = false;
    warningsFound.multNonZEntries = false;
    
    if ~isempty(Uhard(Uhard<0))
        warning('Negative entries found in the hard clustering!\n')
        warningsFound.neg = true;
    end
    if any(all(~Uhard,1))
        warning('Some cluster do not contain any datapoint!\n')
        warningsFound.noDat = true;
    end
    if ~isempty(numEqualAssign(numEqualAssign>1))
        warning('Some datapoints in the hard clustering have the same probability assigned to different centroids.\n')
        warningsFound.prob = true;
    end
    if any(sum(~~Uhard,2)>1)
        warning('There are multiple nonzero entries in the rows of the hard clustering matrix U.')    
        warningsFound.multNonZEntries = true;
    end
%     else
%         cprintf('green', 'No warnings found for the hard clustering!\n')
%     end
end