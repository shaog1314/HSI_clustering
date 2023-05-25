clear all
nC=10;
lam=1e3;

load(['data.mat']);% X is the 3-D HSI cube and grd is the ground truth
l=unique(grd);
lnum=l(l~=0);
for i=1:length(lnum)
    grd(grd==lnum(i))=i;
end

[row,col,band]=size(X);
X_vec=reshape(double(X),row*col,band)';
Y=data_pre_processing_all(X,grd);


r=0;
rho=1;

%% superpixels segmentation  %%%%%%%%
W =pca(X_vec','NumComponents',1);
X_p=X_vec'*W;
X_pt=reshape(X_p,row,col,1);
X_pt=X_pt-min(X_pt(:));
X_pt=X_pt/max(X_pt(:))*255;

tic
Sup_map1 = mex_ers(X_pt,nC); % segmentation map with matrix of 0 - nC
[out.kappa,out.AA,out.OA,out.ACC,out.NMI,out.ARI,out.PUR,~,out.grps] = Sup_SSC(Y,r,lam,Sup_map1,rho,grd,1);
out.time=toc;