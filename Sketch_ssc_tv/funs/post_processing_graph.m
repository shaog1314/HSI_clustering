function [idx] = post_processing_graph(Z,K,k,sigma,large_scale,step,use_ANN,ANN_comp_prcnt,ANN_num_trees,spectral_method,vl_feat)

if nargin < 4
    sigma = 0;
end
if nargin < 5
    large_scale = true;
end
if nargin < 6
    step = 0;
end

if nargin < 7
    use_ANN = 1;
end
if nargin < 8
    ANN_comp_prcnt = 0.1;
end
if nargin < 9
    ANN_num_trees = 10;
end
if nargin < 10
    spectral_method = 1;
end
if nargin < 11
    vl_feat = 0;
end
if sigma == 0
    gwp = 'default';
else
    gwp = sigma;
end
options = struct('Kernel','linear', ...
                 'KernelParam',1, ...
                 'NN',k,...
                 'k' , 2, ...
                 'GraphDistanceFunction','euclidean',... 
                 'GraphWeights', 'binary', ...
                 'GraphWeightParam',gwp, ...
                 'GraphNormalize',1, ...
                 'ClassEdges',0,...
                 'LaplacianDegree',1, ...
                  'gamma_A',1.0,... 
                 'gamma_I',1.0, ...
                  'mu',0.5, ...
                  'step',step, ...
                  'comp_prcnt', ANN_comp_prcnt, ...
                  'ANN_trees', ANN_num_trees); 
options.GraphWeights = 'heat';
options.GraphDistanceFunction = 'Fast_SqDist';
% options.sigma=t_sigma;
if isa(Z,'single')
    Z = double(Z);
end

if ~use_ANN
    A = create_knn_graph(Z,options);
%     [A] = create_ann_graph_grd(Z');
else
    A = create_ann_graph(Z',options);
end

% disp('running spectral clustering');
if spectral_method == 1
    idx = SpectralClustering(A,K,large_scale,0,vl_feat);
elseif spectral_method == 2
    idx = SpectralClustering_alt(A,K,large_scale,0,vl_feat);
end



end

