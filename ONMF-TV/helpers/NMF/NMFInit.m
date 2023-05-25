function [U,V,warningsUInit, warningsVInit] = NMFInit(Y,p,n,m,epsStabInit,iniType)
%Computes the initialization for the NMF algorithm.

%INPUT:
%Y:                 Data matrix
%p:                 Number of extracted Features
%n:                 Number of rows of the data matrix
%m:                 Number of columns of the data matrix
%espStabInit:       Projection value to avoid numerical instabilities
%iniType:           Initialization type. Can be 'manual', 'random',
%                   'svdbased' (SVD), 'kmeansPP' (k-means++).

%OUTPUT:
%U:                 Cluster membership matrix
%V:                 Centroid matrix
%warningsUInit:     Warnings for the initialization of matrix U
%warningsVInit:     Warnings for the initialization of matrix V


% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    switch iniType

        case 'manual'
              %(TODO)

        case 'random'
              % fprintf('Determining a random initialization.\n');
    %           rng(42);
    %           K0 = abs(randn(n,p))+eps;
    %           X0 = abs(randn(p,m))+eps;
              U = abs(randn(n,p));
              V = abs(randn(p,m));

              projectNegativeAndHighValues(U, epsStabInit);
              projectNegativeAndHighValues(V, epsStabInit);


              % fprintf('Trying to optimize a scaling factor k s. t. k*K0*k*X0~Y:\n');
              % There is an analytical solution for the optimal scaling factor:
    %           k_analytical=sqrt(sum(sum(Y.*(K0*X0)))/sum(sum((K0*X0).^2)));
    %           if k_analytical<=0 || k_analytical>1
    %               fprintf('Warning: k_analytical not in [0,1].\n');
    %           end
    %           fprintf('Calculated factor: k_analytical=%6.4f\n',k_analytical); 
    %           K0=k_analytical*K0;
    %           X0=k_analytical*X0;
    %           B=K0;
    %           C=X0;


        case 'svdbased'

    %           fprintf('Calculating an svd based initialization.\n');
    %           fprintf('Taking use of a fast svd variant.\n');
%             rng(42, 'twister');
            %tic; 
            [U,S,V] = MSFsvd(Y,2*p,1,0); % Randomization errors become larger for the 
            %           [U,S,V] = MSFsvd(Y,p+1,1,0); % Randomization errors become larger for the 
            % last calculated svd-vectors depending on the dimension parameter p, 
            % a general recommendation is to calculate a fast svd for
            % roughly 1.5*p vectors and only take the first p calculated vectors.
            % Here, we choose 2*p instead of p.
            %toc;     
            K0 = zeros(n,p);
            X0 = zeros(p,m);
            K0(:,1) = sqrt(S(1,1))*abs(U(:,1));         
            X0(1,:) = sqrt(S(1,1))*abs(V(:,1)'); 
            % For the following compare the publication of Boutsidis/Gallopoulos:
            % SVD based initialization: A head start for nonnegative matrix
            % factorization, in Pattern Recognition 41 (2008), 1350-1362
            posf = @(X)(X>=0).*X;
            negf = @(X)(X<0).*(-X);
            for j=2:p
            uu = U(:,j); vv = V(:,j);
            uup = posf(uu); uun = negf(uu) ;
            vvp = posf(vv); vvn = negf(vv);
            n_uup = norm(uup);
            n_vvp = norm(vvp) ;
            n_uun = norm(uun) ;
            n_vvn = norm(vvn) ;
            termp = n_uup*n_vvp; termn = n_uun*n_vvn;
            if (termp >= termn)
                K0(:,j) = sqrt(S(j,j)*termp)*uup/n_uup; 
                X0(j,:) = sqrt(S(j,j)*termp)*vvp'/n_vvp;
            else
                K0(:,j) = sqrt(S(j,j)*termn)*uun/n_uun; 
                X0(j,:) = sqrt(S(j,j)*termn)*vvn'/n_vvn;
            end
            end



    %         if min(K0(:)) < epsStabInit
    %           K0 = K0+epsStabInit;
    %         end
    %         if min(X0(:)) < epsStabInit
    %           X0 = X0+epsStabInit;

            U=K0;
            V=X0;
%            U = projectNegativeAndHighValues(U, epsStabInit);
%            V = projectNegativeAndHighValues(V, epsStabInit);

        case 'kmeansPP'
            [U, V, ~] = kmeansPP(Y, p);
%            U = projectNegativeAndHighValues(U, epsStabInit);
%            V = projectNegativeAndHighValues(V, epsStabInit);
            
        otherwise
            error('Unknown Initialization Method')
    end
    
    %Check Initialization of the clustering
    [U, V, warningsUInit, warningsVInit]=checkClusteringInit(U, V);
    
    
end