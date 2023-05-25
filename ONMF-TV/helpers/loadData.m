function data = loadData(whichData, dataPath)
%Load the data based on the specified string 'datasetName'.

%INPUT:
%whichData:         Should be either 'exampleData' or 'ownData'. In the
%                   case of 'exampleData', the Indian Pines dataset is
%                   used, which is included in the Image Processing
%                   Toolbox™ Hyperspectral Imaging Library of Matlab. In
%                   the case of 'ownData', the own data of the user is
%                   used (see README.md).
%dataPath:          In case of 'ownData', the full path to the *.mat file
%                   with the dataset and all necessary components is
%                   needed (see README.md).

%OUTPUT:
%data:              struct, with the following fields:
%X:                 Datamatrix: X has to be arranged in such a way, that
%                   the observations are ordered rowwise.
%labels:            Ground Truth of the clustering of the dataset including
%                   the zero spectra as an image. Classes are labelled as
%                   1:NC with NC number of classes
%indexGrid:         x and y dimension of the images, which are gathered in
%                   the columns of U
%indexNonzero:      The vector containing the index values of the non-zero
%                   spectra
%imgMask:           Struct. Image mask, which define the non-annotated
%                   regions or zero-spectra of the dataset for later
%                   plotting.

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.
    
    if nargin < 1 || (nargin < 2 && strcmp(whichData, 'ownData') )
        error('Wrong number of input arguments')
    end
    
    if (strcmp(whichData, 'ownData'))
        data = load(dataPath, 'X', 'labels', 'indexGrid', 'indexNonzero',...
            'imgMask');
    elseif (strcmp(whichData, 'exampleData'))
        labels = load('indian_pines_gt.mat');
        %Get the labels (see README.md)
        data.labels = labels.indian_pines_gt;
        %Get the indexGrid (see README.md)
        data.indexGrid = data.labels>0;
        %Get the indexNonzero Vector (see README.md)
        data.indexNonzero = getIndexNonzeroFromIndexGrid(data.indexGrid);
        %Get the Data as a hypercube 
        hcube = hypercube('indian_pines.dat');
        %Get some optional metaData
        data.metaData = hcube.Metadata;
        %Get the image mask (see README.md)
        data.imgMask.labels = data.indexGrid;
        data.imgMask.zeroSp = ones(size(hcube.DataCube, 1),...
            size(hcube.DataCube, 2));
        %Get the data matrix (see README.md)
        X = reshape(hcube.DataCube, [size(hcube.DataCube, 1)*...
            size(hcube.DataCube, 2), size(hcube.DataCube, 3)]);
        %Restrict the data to the annotated spectra (see README.md)
        X = X(data.indexNonzero, :);
        data.X = X;
    else
        error('The first argument should be either ''exampleData'' or ''ownData'' ')
    end
end
