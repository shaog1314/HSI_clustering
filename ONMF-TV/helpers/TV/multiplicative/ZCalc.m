function [Z_column_Zerospectra]=ZCalc(K_column_reshaped_Zerospectra, nRows, nColumns, epsilon_TV, P_column_Zerospectra)
%Calculate the columns of the matrix Z. For details, see also P. Fernsel,
%P. Maass, A Survey on Surrogate Approaches to Non-negative Matrix
%Factorization, Vietnam Journal of Mathematics, 46, 987-1021, 2018.

%K_column_reshaped_Zerospectra:    Column of the matrix K (with zerospectra)
%nRows:                                                    height of the images
%nColumns:                                               width of the images
%epsilon_TV:                                            parameter of the TV-penalty
%P_column_Zerospectra                       The current column of the matrix P (with zerospectra)

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

Z_column_reshaped_Zerospectra=zeros(nRows, nColumns);
P_column_reshaped_Zerospectra=reshape(P_column_Zerospectra, [nRows nColumns]);

K=K_column_reshaped_Zerospectra;

%Calculation of the boundary values

    %Matrix on the corners

        %Top left corner
        Result_1=1./(sqrt((K(1,1)-K(1,2)).^2 + (K(1,1)-K(2,1)).^2 + epsilon_TV^2));
        Result_2=(K(1,1)+K(1,2))./2 + (K(1,1)+K(2,1))./2;
        Z_column_reshaped_Zerospectra(1,1)=(1./P_column_reshaped_Zerospectra(1,1)).*(Result_1.*Result_2);
        %Bottom left corner
        Result_1=1./(sqrt((K(nRows, 1)-K(nRows, 2)).^2 + epsilon_TV^2));
        Result_2=(K(nRows, 1)+K(nRows, 2))./2;
        Result_3=(K(nRows, 1)+K(nRows-1, 1))./(2.*(sqrt((K(nRows-1, 1)-K(nRows, 1)).^2 + (K(nRows-1, 1)-K(nRows-1, 2)).^2 + epsilon_TV^2)));
        Z_column_reshaped_Zerospectra(nRows,1)=(1./P_column_reshaped_Zerospectra(nRows,1)).*((Result_1.*Result_2) + Result_3);
        %Top right corner
        Result_1=1./(sqrt((K(1, nColumns)-K(2, nColumns)).^2 + epsilon_TV^2));
        Result_2=(K(1, nColumns)+K(2, nColumns))./2;
        Result_3=(K(1, nColumns)+K(1, nColumns-1))./(2.*(sqrt((K(1, nColumns-1)-K(1, nColumns)).^2 + (K(1, nColumns-1)-K(2, nColumns-1)).^2 + epsilon_TV^2)));
        Z_column_reshaped_Zerospectra(1, nColumns)=(1./P_column_reshaped_Zerospectra(1, nColumns)).*((Result_1.*Result_2) + Result_3);
        %Bottom right corner
        Result_1=(K(nRows, nColumns)+K(nRows, nColumns-1))./(2.*(sqrt((K(nRows, nColumns-1)-K(nRows, nColumns)).^2 + epsilon_TV^2)));
        Result_2=(K(nRows, nColumns)+K(nRows-1, nColumns))./(2.*(sqrt((K(nRows-1, nColumns)-K(nRows, nColumns)).^2 + epsilon_TV^2)));
        Z_column_reshaped_Zerospectra(nRows,nColumns)=(1./P_column_reshaped_Zerospectra(nRows, nColumns)).*(Result_1 + Result_2);
    
    %Matrix on the sides (without the corners)
        
        %Left side
        Result_1=1./(sqrt((K(2:nRows-1, 1)-K(3:nRows, 1)).^2 + (K(2:nRows-1, 1)-K(2:nRows-1, 2)).^2 + epsilon_TV^2));
        Result_2=(K(2:nRows-1, 1)+K(3:nRows, 1))./2 + (K(2:nRows-1, 1)+K(2:nRows-1, 2))./2;
        Result_3=(K(2:nRows-1, 1)+K(1:nRows-2, 1))./(2.*(sqrt((K(1:nRows-2, 1)-K(2:nRows-1, 1)).^2 + (K(1:nRows-2, 1)-K(1:nRows-2, 2)).^2 + epsilon_TV^2)));
        Z_column_reshaped_Zerospectra(2:nRows-1,1)=(1./P_column_reshaped_Zerospectra(2:nRows-1, 1)).*((Result_1.*Result_2) + Result_3);
        %Top side
        Result_1=1./(sqrt((K(1, 2:nColumns-1)-K(1, 3:nColumns)).^2 + (K(1, 2:nColumns-1)-K(2, 2:nColumns-1)).^2 + epsilon_TV^2));
        Result_2=(K(1, 2:nColumns-1)+K(1, 3:nColumns))./2 + (K(1, 2:nColumns-1)+K(2, 2:nColumns-1))./2;
        Result_3=(K(1, 2:nColumns-1)+K(1, 1:nColumns-2))./(2.*(sqrt((K(1, 1:nColumns-2)-K(1, 2:nColumns-1)).^2 + (K(1, 1:nColumns-2)-K(2, 1:nColumns-2)).^2 + epsilon_TV^2)));
        Z_column_reshaped_Zerospectra(1, 2:nColumns-1)=(1./P_column_reshaped_Zerospectra(1, 2:nColumns-1)).*((Result_1.*Result_2) + Result_3);
        %Right side
        Result_1=1./(sqrt((K(2:nRows-1, nColumns)-K(3:nRows, nColumns)).^2 + epsilon_TV^2));
        Result_2=(K(2:nRows-1, nColumns)+K(3:nRows, nColumns))./2;
        Result_3=(K(2:nRows-1, nColumns)+K(1:nRows-2, nColumns))./(2.*(sqrt((K(1:nRows-2, nColumns)-K(2:nRows-1, nColumns)).^2 + epsilon_TV^2)));
        Result_4=(K(2:nRows-1, nColumns)+K(2:nRows-1, nColumns-1))./(2.*(sqrt((K(2:nRows-1, nColumns-1)-K(2:nRows-1, nColumns)).^2 + (K(2:nRows-1, nColumns-1)-K(3:nRows, nColumns-1)).^2 + epsilon_TV^2)));
        Z_column_reshaped_Zerospectra(2:nRows-1, nColumns)=(1./P_column_reshaped_Zerospectra(2:nRows-1, nColumns)).*((Result_1.*Result_2) + Result_3 + Result_4);
        %Bottom side
        Result_1=1./(sqrt((K(nRows, 2:nColumns-1)-K(nRows, 3:nColumns)).^2 + epsilon_TV^2));
        Result_2=(K(nRows, 2:nColumns-1)+K(nRows, 3:nColumns))./2;
        Result_3=(K(nRows, 2:nColumns-1)+K(nRows, 1:nColumns-2))./(2.*(sqrt((K(nRows, 1:nColumns-2)-K(nRows, 2:nColumns-1)).^2 + epsilon_TV^2)));
        Result_4=(K(nRows, 2:nColumns-1)+K(nRows-1, 2:nColumns-1))./(2.*(sqrt((K(nRows-1, 2:nColumns-1)-K(nRows, 2:nColumns-1)).^2 + (K(nRows-1, 2:nColumns-1)-K(nRows-1, 3:nColumns)).^2 + epsilon_TV^2)));
        Z_column_reshaped_Zerospectra(nRows, 2:nColumns-1)=(1./P_column_reshaped_Zerospectra(nRows, 2:nColumns-1)).*((Result_1.*Result_2) + Result_3 + Result_4);

        
    %Finally the calculation of the interior values. So we don't have to
    %consider any special cases for the boundary values.

        %Calculation of the "normal" neighborhood
        K_central=K(2:(nRows-1), 2:(nColumns-1));
        K_right=K(2:(nRows-1), 3:nColumns);
        K_bottom=K(3:nRows, 2:(nColumns-1));

        Result_1=1./(sqrt((K_central-K_right).^2 + (K_central-K_bottom).^2 + epsilon_TV^2));
        Result_2=(K_central + K_right)./2 + (K_central  + K_bottom)./2;

        %Calculation of the "adjoint" neighborhood
        K_top=K(1:(nRows-2), 2:(nColumns-1));
        K_left=K(2:(nRows-1), 1:(nColumns-2));
        K_left_bottom=K(3:nRows, 1:(nColumns-2));
        K_top_right=K(1:(nRows-2), 3:nColumns);


        Result_3=(K_central+K_left)./(2.*(sqrt((K_left-K_central).^2 + (K_left-K_left_bottom).^2 + epsilon_TV^2)));
        Result_4=(K_central+K_top)./(2.*(sqrt((K_top-K_central).^2 + (K_top-K_top_right).^2 + epsilon_TV^2)));

    %Final Result
        Z_column_reshaped_Zerospectra(2:(nRows-1), 2:(nColumns-1))=(1./P_column_reshaped_Zerospectra(2:(nRows-1), 2:(nColumns-1))).*((Result_1.*Result_2)+Result_3+Result_4);

    %Reshape the matrix to a column vector
        Z_column_Zerospectra=reshape(Z_column_reshaped_Zerospectra, [nRows*nColumns, 1]);
