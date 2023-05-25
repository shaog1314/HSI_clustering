% Before running the code, VLFeat package is required. see https://www.vlfeat.org/
% Note that after installing VLFeat package, all "mex" related paths under
% vl_toolbox should be moved to top. Otherwise, the k-means function in
% vltoolbox cannot be used.

clear
load('Indian_part.mat')
%% prameters to be specified
lam=1e-3; % penalty of sparsity constraint
lam_tv=1e-2; % penalty of spatial constraint
ro=1; % parameter in the optimization algorithm
n=70; % dictionary size
k_nn=30; % the number of nearest neighbors in knn/ann graph
%% other parameters
Rand_mat_type = 'Rademacher'; %type of random matrix to use;
use_ANN = 1; %Set to 1 to use ANN methods (requires VLFeat package)
ANN_num_trees = ceil(n/21); %Parameters for ANN computations
ANN_comp_prcnt = 0.05;
spectral_method = 1;
opts.lam=lam;
opts.lam_tv=lam_tv;
opts.tol=1e-5;
opts.max_itr=100;
opts.ro=ro;
opts.step=1.1;
vl_feat = true; %Set to true if you want to use the VLfeat package 
large_scale_spectral = true;
knn_graph_step = 500; 


[row,col,band]=size(X);
opts.row=row;
opts.col=col;
opts.band=band;
N=row*col;
X_vec=reshape(X,row*col,band)';
[Y]=normc(X_vec);   % normalize input data (each column is a data point) in l2 norm
No_classes=max(grd(:)); % number of clusters
%% obtain random projection matrix
R = construct_random_mat(N,n,Rand_mat_type); %Generate random matrix
D = Y*R; %form sketched dictionary
%% get coefficient matrix
Z=sketch_tv_sc(Y,D,opts);
Z = normc(Z);
%% Generation of k-nn matrix and spectral clustering.
grps = post_processing_graph(Z',No_classes,k_nn,0,large_scale_spectral,knn_graph_step,use_ANN,ANN_comp_prcnt,ANN_num_trees,spectral_method,vl_feat);
%% Evaluation of clustering
[kappa,AA,OA,ACC,grp] = fast_performance(grd(:),grps);


                           
                           
                           
                           
                           
                          