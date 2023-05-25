function Y_t=data_pre_processing_all(X,grd)
% extract data with spatial window after applying pca
M_p=3;
M_z=2*M_p+1;
[row,col,band]=size(X);
X_vec=reshape(double(X),row*col,band)';
pcs=8;
W =pca(X_vec','NumComponents',pcs);
X_p=X_vec'*W;
X_pt=reshape(X_p,row,col,pcs);

img_2=padarray(X_pt,[M_p M_p 0],'symmetric'); %normalize or not
[ind_row, ind_col]=find(grd~=-1);
ind_row_1=ind_row+M_p;
ind_col_1=ind_col+M_p;
Y=zeros(length(ind_row_1),pcs*M_z^2);
for i=1:length(ind_row_1)
    Y(i,:)=reshape(img_2(ind_row_1(i)-M_p:ind_row_1(i)+M_p,ind_col_1(i)-M_p:ind_col_1(i)+M_p,:),1,pcs*M_z^2);
end
Y_t=normr(Y)';
% Y_t=reshape(normr(Y),row,col,size(Y,2));