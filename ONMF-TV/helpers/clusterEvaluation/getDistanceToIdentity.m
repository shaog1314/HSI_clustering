function [distUtoIl1, distUtoIl2, distUtoImax]=getDistanceToIdentity(U)
%Compute the l1, l2 distance and the maximum absolute value of the matrix
%U-I of the cluster membership matrix U to the identity matrix I.

%INPUT:
%U:                 Cluster Membership Matrix
%measure:           String. Can be 'l1' or 'l2' for l1 distance measure and
%                   l2 distance measure respectively.

%OUTPUT:
%distUtoIl1:          l1 distance
%distUtoIl2:          l2 distance
%distUtoImax:         Maximum value


% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    distUtoIl1=sum(vecnorm(U'*U-eye(size(U,2)),1,1));
    distUtoIl2=norm(U'*U-eye(size(U,2)), 'fro');
    distUtoImax=max(max(abs(U'*U-eye(size(U,2)))));
    
end