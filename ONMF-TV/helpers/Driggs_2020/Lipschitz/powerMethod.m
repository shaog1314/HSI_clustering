function [lambda, EV]=powerMethod(U, V, W, param, whichVar)
%Estimation of the Lipschitz constant of the gradient of the NMF cost
%function (based on Driggs, D.; Tang, J.; Liang, J.; Davies, M.; Schönlieb,
%C.B. SPRING: A fast stochastic proximal alternating method for non-smooth
%non-convex optimization. arXiv preprint 2020, arXiv: 2002.12266) for PALM.
%We consider the NMF problem
%1/2*||X-U*V||_F^2 + sigma1/2 * ||W^T*U - I||_F^2 + sigma2/2*||W-U||_F^2.

%INPUT:
%U:                 Cluster membership matrix
%V:                 Matrix, which has the centroids in its rows
%W:                 Auxiliary Variable of the size of U
%param:             Struct with all hyperparameters
%whichVar:          String, which specifies the variable of the
%                   corresponding partial derivative. U for df/dU, V for
%                   dF/dV and W analogously, if F is the NMF cost function.

%OUTPUT:
%lambda:            Maximal Eigenvalue (and therefore the estimated
%                   Lipschitz constant of the gradient.
%                   In the case of
%                   whichVar='U', lambda is a vector of two EV [lambda1,
%                   lambds 2], where lambda1 is the EV of V*V' and lambda 2
%                   of param.sigma1*W*W' + param.sigma2*eye(M). The
%                   corresponding eigenvectors are stored columnwise in EV.
%EV:                Optional: A corresponding eigenvector(s) of the
%                   eigenvalue(s) (see description of lambda).

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    if nargin ~= 5
        error('Wrong number of input arguments')
    end
    
    %Define the dimensions of the matrices according to the notation
    K = size(U,2);
%     N = size(V, 2);
    
    switch whichVar
        case 'U' %In this case, we have to run the power method two times.
            %For V*V'
            xi = randn(K, 1);
            xi = xi./norm(xi);
            VV = V*V';
            for i=1:param.powerIt
                x_tmp = VV*xi;
                xi = x_tmp./norm(x_tmp);
            end
            EV1 = xi;
            lambda1 = norm(VV * EV1);
            %For param.sigma1*W*W' + param.sigma2*eye(M)
            xi = randn(K, 1);
            xi = xi./norm(xi);
            WW = param.sigma1 * (W'*W) + param.sigma2*eye(K);
            for i=1:param.powerIt
                x_tmp = WW*xi;
                xi = x_tmp./norm(x_tmp);
            end
            EV2 = xi;
            lambda2 = norm(WW * EV2);
            %Save both eigenvalues and eigenvectors to the output variables
            lambda = [lambda1, lambda2];
            EV = [EV1, EV2];
        case 'V'
            xi = randn(K, 1);
            xi = xi./norm(xi);
            UU = U'*U;
            for i=1:param.powerIt
                x_tmp = UU*xi;
                xi = x_tmp./norm(x_tmp);
            end
            EV = xi;
            lambda = norm(UU * EV);
        case 'W'
            xi = randn(K, 1);
            xi = xi./norm(xi);
            UU = param.sigma1 * (U'*U) + param.sigma2*eye(K);
            for i=1:param.powerIt
                x_tmp = UU*xi;
                xi = x_tmp./norm(x_tmp);
            end
            EV = xi;
            lambda = norm(UU * EV);
        otherwise
            error('Unknown variable of corresponding partial derivative')
    end

end