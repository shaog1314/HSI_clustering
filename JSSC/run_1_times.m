function run_1_times(l1,l2,nCs,lams,Sup_map1,Y,r,rho,grd,mu,savename)
for i=1:l1
    nC=nCs(i);
    tic
%     Sup_map = mex_ers(X_pt,nC); % segmentation map with matrix of 0 - nC
    Sup_map=Sup_map1(:,:,i);
    %% segmentation show
% [height,width] = size(X_pt);
% [bmap] = seg2bmap(Sup_map,width,height);
% img=X_pt;
% grey_img=img;
% bmapOnImg = img;
% idx = find(bmap>0);
% timg = grey_img;
% timg(idx) = 255;
% pimg=grey_img;
% pimg(idx)=0;
% qimg=grey_img;
% qimg(idx)=0;
% bmapOnImg(:,:,2) = timg;
% bmapOnImg(:,:,1) = pimg;
% bmapOnImg(:,:,3) = qimg;
% bmapOnImg=bmapOnImg/255;
% figure;subplot(121); imagesc(grd);axis image;
% subplot(122); imshow(bmapOnImg,[]);
%%

    for j=1:l2
        lambda=lams(j);
        [kappa(i,j),AA(i,j),OA(i,j),ACC(:,i,j),~,grps(:,i,j),index] = Sup_SSC(Y,r,lambda,Sup_map,rho,grd,mu);
        max(OA,[],2)
        time(i,j)=toc;
        save(savename,'kappa','AA','OA','ACC','grps','lams','nCs','time');
    end


end