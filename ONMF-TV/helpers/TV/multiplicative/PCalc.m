function [P_column_Zerospectra]=PCalc(K_column_reshaped_Zerospectra, nRows, nColumns, epsilon_TV)
%Calculate the columns of the matrix P. For details, see also P. Fernsel,
%P. Maass, A Survey on Surrogate Approaches to Non-negative Matrix
%Factorization, Vietnam Journal of Mathematics, 46, 987-1021, 2018.

%K_column_reshaped_Zerospectra:       Column of the matrix K (with zerospectra)
%nRows:                               height of the images
%nColumns:                            width of the images
%epsilon_TV:                          parameter of the TV-penalty

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

P_column_reshaped_withZeros=zeros(nRows, nColumns);

K=K_column_reshaped_Zerospectra;

%Calculation of the boundary values

    %Matrix on the corners

        %Top left corner
        P_column_reshaped_withZeros(1,1)=2./(sqrt((K(1,1)-K(1,2)).^2 + (K(1,1)-K(2,1)).^2 + epsilon_TV^2));
        %Bottom left corner
        Result_1=1./(sqrt((K(nRows, 1)-K(nRows, 2)).^2 + epsilon_TV^2));
        Result_2=1./(sqrt((K(nRows-1, 1)-K(nRows, 1)).^2 + (K(nRows-1, 1)-K(nRows-1, 2)).^2 + epsilon_TV^2));
        P_column_reshaped_withZeros(nRows,1)=Result_1 + Result_2;
        %Top right corner
        Result_1=1./(sqrt((K(1, nColumns)-K(2, nColumns)).^2 + epsilon_TV^2));
        Result_2=1./(sqrt((K(1, nColumns-1)-K(1, nColumns)).^2 + (K(1, nColumns-1)-K(2, nColumns-1)).^2 + epsilon_TV^2));
        P_column_reshaped_withZeros(1,nColumns)=Result_1 + Result_2;
        %Bottom right corner
        Result_1=1./(sqrt((K(nRows, nColumns-1)-K(nRows, nColumns)).^2 + epsilon_TV^2));
        Result_2=1./(sqrt((K(nRows-1, nColumns)-K(nRows, nColumns)).^2 + epsilon_TV^2));
        P_column_reshaped_withZeros(nRows,nColumns)=Result_1 + Result_2;
    
    %Matrix on the sides (without the corners)
        
        %Left side
        Result_1=2./(sqrt((K(2:nRows-1, 1)-K(3:nRows, 1)).^2 + (K(2:nRows-1, 1)-K(2:nRows-1, 2)).^2 + epsilon_TV^2));
        Result_2=1./(sqrt((K(1:nRows-2, 1)-K(2:nRows-1, 1)).^2 + (K(1:nRows-2, 1)-K(1:nRows-2, 2)).^2 + epsilon_TV^2));
        P_column_reshaped_withZeros(2:nRows-1,1)=Result_1 + Result_2;
        %Top side
        Result_1=2./(sqrt((K(1, 2:nColumns-1)-K(2, 2:nColumns-1)).^2 + (K(1, 2:nColumns-1)-K(1, 3:nColumns)).^2 + epsilon_TV^2));
        Result_2=1./(sqrt((K(1, 1:nColumns-2)-K(1, 2:nColumns-1)).^2 + (K(1, 1:nColumns-2)-K(2, 1:nColumns-2)).^2 + epsilon_TV^2));
        P_column_reshaped_withZeros(1, 2:nColumns-1)=Result_1 + Result_2;
        %Right side
        Result_1=1./(sqrt((K(2:nRows-1, nColumns)-K(3:nRows, nColumns)).^2 + epsilon_TV^2));
        Result_2=1./(sqrt((K(2:nRows-1, nColumns-1)-K(2:nRows-1, nColumns)).^2 + (K(2:nRows-1, nColumns-1)-K(3:nRows, nColumns-1)).^2 + epsilon_TV^2));
        Result_3=1./(sqrt((K(1:nRows-2, nColumns)-K(2:nRows-1, nColumns)).^2 + epsilon_TV^2));
        P_column_reshaped_withZeros(2:nRows-1, nColumns)=Result_1 + Result_2 + Result_3;
        %Bottom side
        Result_1=1./(sqrt((K(nRows, 2:nColumns-1)-K(nRows, 3:nColumns)).^2 + epsilon_TV^2));
        Result_2=1./(sqrt((K(nRows-1, 2:nColumns-1)-K(nRows, 2:nColumns-1)).^2 + (K(nRows-1, 2:nColumns-1)-K(nRows-1, 3:nColumns)).^2 + epsilon_TV^2));
        Result_3=1./(sqrt((K(nRows, 1:nColumns-2)-K(nRows, 2:nColumns-1)).^2 + epsilon_TV^2));
        P_column_reshaped_withZeros(nRows, 2:nColumns-1)=Result_1 + Result_2 + Result_3;

    %Finally the calculation of the interior values. So we don't have to
    %consider any special cases for the boundary values.

    %Calculation of the "normal" neighborhood
    K_central=K(2:(nRows-1), 2:(nColumns-1));
    K_right=K(2:(nRows-1), 3:nColumns);
    K_bottom=K(3:nRows, 2:(nColumns-1));

    Result_1=2./(sqrt((K_central-K_right).^2 + (K_central-K_bottom).^2 + epsilon_TV^2)); %Erste Teilsumme

    %Calculation of the "adjoint" neighborhood
    K_top=K(1:(nRows-2), 2:(nColumns-1));
    K_left=K(2:(nRows-1), 1:(nColumns-2));
    K_left_bottom=K(3:nRows, 1:(nColumns-2));
    K_top_right=K(1:(nRows-2), 3:nColumns);

    Result_2=1./(sqrt((K_left-K_central).^2 + (K_left - K_left_bottom).^2 + epsilon_TV^2));
    Result_3=1./(sqrt((K_top-K_central).^2 + (K_top - K_top_right).^2 + epsilon_TV^2));

    %Final Result
    P_column_reshaped_withZeros(2:(nRows-1), 2:(nColumns-1))=Result_1 + Result_2 + Result_3;

    %Reshape the matrix to a column vector
    P_column_Zerospectra=reshape(P_column_reshaped_withZeros, [nRows*nColumns, 1]);