function [D] = Fast_SqDist(X,Y)
%FAST_SQDIST Function that computes the euclidean distance between the rows
%of matrices X and Y fast.

D = bsxfun(@plus, dot(X,X,2), dot(Y,Y,2)') - 2 * (X * Y');

D(D <0) = 0;
end

