function [U, V, index] = kmeansPP(X, K)
%Performs the initialization based on the k-means++ algorithm, which uses
%an heuristic to find centroid seeds for k-means clustering. According to
%Arthur and Vassilvitskii [1], k-means++ improves the running time of
%Lloyd’s algorithm, and the quality of the final solution.

%INPUT:
%X:                 Datamatrix of the size M-by-N with M observations and N
%                   features. 
%K:                 Number of Clusters

%OUTPUT:
%U:                 Cluster Membership Matrix of the size M-by-K
%V:                 Matrix of the size K-by-N containing the centroids
%                   row-wise.
%index:             Vector with the indices of the chosen datavectors as
%                   centroids


% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

%     %Check Input
%     if nargin<=2
%         epsStabInit=[1e-4, 10^35];
%     end

    %Choice of the centroids

    S = RandStream.getGlobalStream;
    [M, N] = size(X);
    X = X';
    
    % Select the first seed by sampling uniformly at random
    index = zeros(1,K);
    
    V = zeros(N, K); %Transposed
%     U = zeros(M, K);
    
    [V(:,1), index(1)] = datasample(S,X,1,2);
    minDist = inf(M,1);
    
    % Select the rest of the seeds by a probabilistic model
    for ii = 2:K                    
        minDist = min(minDist,distfun(X,V(:,ii-1)));
        denominator = sum(minDist);
        if denominator==0 || isinf(denominator) || isnan(denominator)
            V(:,ii:k) = datasample(S,X,K-ii+1,2,'Replace',false);
            break;
        end
        sampleProbability = minDist/denominator;
        [V(:,ii), index(ii)] = datasample(S,X,1,2,'Replace',false,...
            'Weights',sampleProbability);        
    end
    
    %Compute the corresponding cluster membership matrix
    D = distfun(X, V);
    [~, idx] = min(D, [], 2);
    U=indexVectorToClusterMembership(idx);
    
    %Transpose V to its original shape
    V=V';
    
    %Set all the zero values in U to 0.5.
%     U(U==0)=0.5; %0.1
    
    U(U==0) = rand(length(U(U==0)),1)*0.5;
    
    %Test for zero entries in the matrices U and V and perform, if needed,
    %a projection step
%     if any(any(U==0))
%         warning('The Matrix U contains zero elements. A suitable projection step will be performed to avoid zero entries.\n')
%         projectNegativeAndHighValues(U, epsStabInit);
%     end
% 
%     if any(all(~U,1))
%         warning('At Initialization, some cluster do not contain any datapoint!\n')
%     elseif ~isempty(numEqualAssign(numEqualAssign>1))
%         warning('At initialization, some datapoints in the hard clustering have the same probability aignment to different centroids.\n')
%     else
%         cprintf('green', 'No warnings found for the initialization of the hard clustering!\n')
%     end
%         
% 
%     if any(any(V==0))
%         warning('The Matrix V contains zero elements. A suitable projection step will be performed to avoid zero entries.\n')
%         projectNegativeAndHighValues(V, epsStabInit);
%     end

end %main function


function D = distfun(X, V)
%Calculates point to cluster centroid distances based on the squared
%euclidean norm.

%INPUT:
%X:                 Datamatrix of the size M-by-N with M observations and N
%                   features. 
%V:                 Matrix of the size K-by-N containing the centroids
%                   row-wise.

%OUTPUT:
%D:                 Value of distance

        D = internal.stats.pdist2mex(X,V,'sqe',[],[],[],[]);
        
end % function