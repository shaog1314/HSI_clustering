function out=run_jssc(data_name,data_type,nC,lam)
switch data_type
    case 'part'
        load(['.\data\' data_name '\data.mat']);
        l=unique(grd);
        lnum=l(l~=0);
        for i=1:length(lnum)
            grd(grd==lnum(i))=i;
        end
    case 'whole'
        load(['.\data\' data_name '\data_whole.mat']);
end
% X=Indian_cut;
% grd=grd_cut;

[row,col,band]=size(X);
X_vec=reshape(double(X),row*col,band)';
Y=data_pre_processing_all(X,grd);


r=0;                     % don't use pca
rho=1;

%% superpixels segmentation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% obtain superpixels
W =pca(X_vec','NumComponents',1);
X_p=X_vec'*W;
X_pt=reshape(X_p,row,col,1);
X_pt=X_pt-min(X_pt(:));
X_pt=X_pt/max(X_pt(:))*255;
%%  Parameters to be changed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nCs=[5 10 20 40 80 160 320];

% lams=[1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4];
tic
Sup_map1 = mex_ers(X_pt,nC); % segmentation map with matrix of 0 - nC
[out.kappa,out.AA,out.OA,out.ACC,out.NMI,out.ARI,out.PUR,~,out.grps,~] = Sup_SSC(Y,r,lam,Sup_map1,rho,grd,1);
out.time=toc;