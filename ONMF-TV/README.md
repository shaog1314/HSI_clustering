# Spatially Coherent Clustering Based on Orthogonal Nonnegative Matrix Factorization

Support code for the following article (hereinafter referred to as "the article"):

P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
Nonnegative Matrix Factorization", Journal of Imaging, 2021.

This paragraph contains a high-level description of the package, with a brief overview of its features and limitations.


## Description

The provided MATLAB® code was used for the numerical evaluation of the numerous orthogonal nonnegative matrix factorization (ONMF) models considered in the above mentioned article. Restrictions apply to the availability of the hyperspectral dataset used in this work as it was obtained with the permission a third party. Hence, this code is written for general nonnegative hyperspectral datasets, which have to be provided by the user. However, an example dataset is provided, which is the 'Indian Pines dataset' coming from Image Processing Toolbox™ Hyperspectral Imaging Library of MATLAB®. For more information, we refer the reader to the chapter 'Providing and Loading the Data'. The results in the article were obtained by executing the provided codes with MATLAB® R2021a.

The evaluation of the different ONMF models is done by executing the corresponding sections in the master script (Evaluation_MasterScript.m). It consists of the following main ONMF algorithms.

### Comparative Methods of the Article

* **KMeans_TV.m** (Classical K-means clustering algorithm + TV Postprocessing)
* **ONMF_TV_Choi.m** (Alternating multiplicative update rules + TV Postprocessing)
* **ONMF_TV_Ding.m** (Alternating multiplicative update rules + TV Postprocessing)
* **ONMF_TV_Pompili1.m** (Alternating expectation-maximization algorithm + TV Postprocessing)
* **ONMF_TV_Pompili2.m** (Alternating algorithm based on an augmented Lagrangian approach with strict orthogonality constraints + TV Postprocessing)
* **ONMF_TV_Kimura.m** (ONMF model + TV Postprocessing)
* **ONMF_TV_Li.m** (ONMF model + TV Postprocessing)

### Proposed Methods of the Article

* **ONMFTV_MUL1.m** (Alternating multiplicative update rules with included TV regularization)
* **ONMFTV_MUL2.m** (Alternating multiplicative update rules with included TV regularization)
* **ONMFTV_PALM.m** (Proximal Alternating Linearized Minimization scheme with included TV regularization)
* **ONMFTV_iPALM.m** (Inertial Proximal Alternating Linearized Minimization scheme with included TV regularization)
* **ONMFTV_SPRING.m** (Stochastic Proximal Alternating Linearized Minimization scheme with included TV regularization)

### Providing and Loading the Data

As described above, the data has to be provided by the user. To do so, the full path to the respective *.mat file has to be provided in the master script (see the variables 'whichData' and 'dataPath'). The corresponding *.mat file has to provide the following variables with the following names (see also loadData.m):

* **indexGrid**: Logical x-by-y matrix, where x and y are the dimensions of the images, which are induced of the channels of the hyperspectral dataset (columns of X or U). An entry of this matrix is one, if a spectrum in the respective position was acquired **and** is labelled. If no spectrum was acquired at the respective position ('zero spectra') or if it is not labelled, the respective entry is zero.
* **labels**: x-by-y matrix (same dimension as indexGrid), which provides the ground truth of the clustering of the dataset as an image including the zero spectra. The classes are labelled as 1:NC with NC being the number of classes (see also Figure 1(b) in the article). It is the same matrix as indexGrid, just with the additional information, which position corresponds to which class.
* **X**: M-by-N matrix of the nonnegative hyperspectral dataset. The matrix should be arranged in such a way, that the rows of the matrix correspond to the measured spectra. The channels of the dataset are stored column-wise in the matrix. The rows of the matrix does **not** include the zero spectra.
* **indexNonzero**: M-by-1 vector (same dimension as the number of rows of X). This vector contains the index values of the non zero spectra.
* **imgMask**: Struct with the following two fields:
	* **zeroSp**: Logical x-by-y matrix (same dimension as indexGrid). This matrix may differ slightly from indexGrid: An entry of this matrix is one, if a spectrum in the respective position was acquired independent from its labelling. If no spectrum was acquired at the respective position ('zero spectra'), the respective entry is zero.
	* **labels**: Logical x-by-y matrix (same dimension as indexGrid). The same matrix as indexGrid.
	
However, as described above, an example dataset is also provided within the code. The used data is the 'Indian Pines dataset' coming from Image Processing Toolbox™ Hyperspectral Imaging Library of MATLAB®. The dataset consists of a single hyperspectral image of size 145-by-145 pixels with 220 color channels. The dataset also contains a ground truth label image with 16 classes. Note that a MATLAB® version of R2020a or later is needed to run the code with the example data (see also loadData.m).

### Choice of Hyperparameters

Note that the used regularization parameters in the codes are optimized with respect to the hyperspectral dataset used in the article (unless otherwise stated in the code). To obtain optimal results for the dataset provided by the user, parameter tests of all hyperparameters are necessary. This is especially true for the provided example dataset, which is described above. In this respect, we note that the comparative method ONMF-TV-Li leads to numerical instabilities with the used regularization parameters.
Furthermore, note that the code comes with no guarantee or warranty of any kind.

## Authors and contributors

**Pascal Fernsel** (Center for Industrial Mathematics, University of Bremen)
Email: p.fernsel@uni-bremen.de
