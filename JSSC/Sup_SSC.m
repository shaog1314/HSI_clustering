function [kappa,AA,OA,ACC,NMI,ARI,PUR,C,grps] = Sup_SSC(X,r,lambda,Sup_map,rho,grd,mu)

p=find(grd~=0);
s=grd(p);

n = length(unique(s));
Xp = DataProjection(X,r);
C=admm_funcs_Sup_SSC_stand(Xp,lambda,Sup_map,mu);
%%
CKSym = BuildAdjacency(thrC(C,rho),0,1);
disp('adjaceny completed.....');
grps = SpectralClustering_ssc(CKSym,n);
[kappa, AA,OA,ACC,NMI,ARI,PUR,grps] = fast_performance_overview(grd(:),grps(:));
