%--------------------------------------------------------------------------
% This function takes a NxN coefficient matrix and returns a NxN adjacency
% matrix by choosing the K strongest connections in the similarity graph
% CMat: NxN coefficient matrix
% K: number of strongest edges to keep; if K=0 use all the exiting edges
% CKSym: NxN symmetric adjacency matrix
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------


function [CKSym,CAbs] = BuildAdjacency(CMat,K,sk)

if (nargin < 2)
    K = 0;
end

N = size(CMat,1);
CAbs = abs(CMat);

[Srt,Ind] = sort( CAbs,1,'descend' );

if (K == 0)
    for i = 1:N
        CAbs(:,i) = CAbs(:,i) ./ (CAbs(Ind(1,i),i)+eps);
    end
else
    for i = 1:N
        for j = 1:K
            CAbs(Ind(j,i),i) = CAbs(Ind(j,i),i) ./ (CAbs(Ind(1,i),i)+eps);
        end
    end
end
ctc=CAbs'*CAbs;
cc=CAbs*CAbs;
switch sk
    case 1
        CKSym = CAbs + CAbs';
    case 2
        CKSym = ctc+ctc';
    case 3
        CKSym = cc+cc';
    case 4
        CKSym = ctc;
    case 5
        CKSym=ctc';
    case 6
        CKSym=ctc+CAbs + CAbs';
    case 7
        CKSym=ctc'+CAbs + CAbs';
    case 8
        CKSym = ctc+ctc'+CAbs + CAbs';
    case 9
        CKSym = cc+cc'+CAbs + CAbs';
    case 10
        CKSym = ctc+cc+cc'+CAbs + CAbs';
    case 11
        CKSym = ctc'+cc+cc'+CAbs + CAbs';
    case 12
        CKSym = ctc+ctc'+cc+cc';
    case 13
        CKSym = ctc+ctc'+cc+cc'+CAbs + CAbs';
end