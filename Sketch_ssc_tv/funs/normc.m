function [X_norm] = normc(X)

X_norm = normr(X');
X_norm = X_norm';
end

