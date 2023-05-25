function Uhard=convertToHardClustering(Ufuzzy)
%Converts the (fuzzy) cluster membership matrix U to a matrix with hard
%clusters by using the maximum across the rows of U.

%INPUT:
%Ufuzzy:            Fuzzy Cluster Membership Matrix

%OUTPUT:
%Uhard:             Matrix with hard clusters

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

%     if size(Ufuzzy,1)<=size(Ufuzzy,2)
%         warning('Probably wrong orientation of the Cluster Membership Matrix!')
%     end
    Ufuzzy(Ufuzzy~=repmat(max(Ufuzzy,[],2),1,size(Ufuzzy,2))) = 0;
    
    
    %Consider now the case if there are equal numbers in some rows of
    %Ufuzzy
    numEqualAssign=sum(Ufuzzy==max(Ufuzzy,[],2),2);
%     if ~isempty(numEqualAssign(numEqualAssign>1))
%         warning('Some datapoints have the same probability aignment to different centroids.')
%     end
    
    dataPointIndexEqualProb = find(numEqualAssign>1);
%     [~,maxInd] = max(Ufuzzy,[],2);
%     maxInd(dataPointIndexEqualProb)
    
    posEqualProb=Ufuzzy==max(Ufuzzy,[],2);
    
    
%     posEqualProb=posEqualProb(dataPointIndexEqualProb,:);
    
    %For reproducability, set the seed:
%     rng(39, 'twister');
    chosenPos=zeros(length(dataPointIndexEqualProb),1);
    for i=1:length(dataPointIndexEqualProb)
        posEqualProbRow=find(posEqualProb(dataPointIndexEqualProb(i),:));
        randomPos=randi(length(posEqualProbRow));
        chosenPos(i)=posEqualProbRow(randomPos);
        posEqualProb(dataPointIndexEqualProb(i),:)=0;
        posEqualProb(dataPointIndexEqualProb(i),chosenPos(i))=1;
    end
    
    
    Uhard=posEqualProb;
    
%     if sum(~any(Uhard,1))>0
%         warning('There exist some clusters, which do not have any associated datapoints.')
%     end
    
end