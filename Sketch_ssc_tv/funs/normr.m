function [nA] = normr(A)


nA = bsxfun(@rdivide,A,sqrt(sum(A.^2,2)));
nA(~isfinite(nA)) = 1; 
end

