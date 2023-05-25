function idx = SpectralClustering(CKSym,K,large_scale,LSR,vl_feat)

% warning off;
if nargin <3 
    large_scale = true;
end
if nargin < 4
    LSR = false;
end
N = size(CKSym,1);
MAXiter = 1000; % Maximum number of iterations for KMeans 
REPlic = 20; % Number of replications for KMeans

% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}

DN = sum(CKSym,2);


DN(DN~=0) = (1./sqrt(DN(DN~=0)));
DN=spdiags(DN,0,speye(size(CKSym,1)));
LapN = speye(N) - DN * CKSym * DN;

clear CKSym
if large_scale
    diff = eps;
    try
        [Vn,~] = eigs(LapN,K,'sr');  % first option
        
    catch ME
        ME
        disp('eigs failed, retrying');
        
        LapN = LapN + 10*eps*speye(N);
        [Vn,~] = eigs(LapN,K,diff);
        
    end
else
    [~,~,vN] = svd(full(LapN));
    Vn = vN(:,N-K+1:N);
end

if ~LSR
    U = bsxfun(@rdivide, Vn, sqrt(sum(Vn.^2, 2))+eps);
else
    U = DN*Vn;
end
clear DN Vn LapN

if vl_feat
    [~,idx] = vl_kmeans(U',K,'MaxNumIterations',MAXiter,'NumRepetitions',REPlic);
else
    idx = kmeans(U,K,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
end