function X = projectOnPositiveValues(X, epsLim)
%Projects all values smaller than epsLim in the given matrix to epsLim.

%INPUT:
%X:             Given Matrix
%epsLim:        Limit of the Projection Values

%OUTPUT
%X:             Projected Matrix

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

  if min(X(:))<epsLim(1)
      X(X<epsLim(1)) = epsLim(1);
  end
end