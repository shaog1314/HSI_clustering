function [A] = create_ann_graph(X,options)
%CREATE_ANN_GRAPH Function that creates a k-nearest neighbor graph using
%approximate nearest neighbor methods (requires VLFeat package) 


NN = options.NN;
num_trees = options.ANN_trees;
n = size(X,2);

% fprintf(1, ...
% ['Computing ' num2str(options.NN) '-ANN  Adjacency Graph \n']);


if options.comp_prcnt < 1
    MaxComp = ceil(options.comp_prcnt*n);
else
    MaxComp = ceil(options.comp_prcnt);
end

% disp('Generating KDtree');
kdtree = vl_kdtreebuild(X,'NumTrees',num_trees);

% disp('Querying tree');
[i,dists] = vl_kdtreequery(kdtree,X,X,'NumNeighbors',NN+1,'MaxComparisons',MaxComp); 
% disp('Creating graph');

qq = double(i(:)); 
I=repmat((1:n),[options.NN+1 1]); I=I(:);

switch options.GraphWeights
    case 'binary'
        A=sparse(I,qq,1,n,n);
        
    case 'heat'
        ZZ = dists(:); 
        t = mean(ZZ);
%%         t = 2*t^2; original version
        t=2*t;   % version of May 25
%         t=options.sigma^2;
        ZZ_exp = exp(-ZZ./t);
        A=sparse(I,qq,ZZ_exp,n,n);



%         ZZ = dists(:); 
%         t = mean(sqrt(ZZ));
%         t=10*t;
%         t=2*t^2;    %??????????
%         ZZ_exp = exp(-ZZ./t);
%         A=sparse(I,qq,ZZ_exp,n,n);
    otherwise
    		 error('Unknown Weighttype\n');   
end
A=A+((A~=A').*A');
end

