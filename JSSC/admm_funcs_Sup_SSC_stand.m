function C=admm_funcs_Sup_SSC_stand(Y,lambda,Sup_map,mu)
MaxItr=30; 
tol=1e-7;
u=mu;
ro=1.1;
u_max=1e6;
[m,n]=size(Y);
C=zeros(n,n);
A=zeros(n,n);
Y1=zeros(1,n);
Y2=zeros(n,n);
yty=Y'*Y;
Sup_map=Sup_map+1;
nC=max(Sup_map(:));
for i=1:nC
    wg(i)=sqrt(sum(sum(Sup_map==i)))/sqrt(nC);
end
wg=wg./sum(wg);


for i=1:MaxItr
    C=(lambda*yty+u*(eye(n)+1))\(lambda*yty+u*(A+1)-ones(n,1)*Y1-Y2);
    for j=1:nC
        idx=find(Sup_map==j);
        A(:,idx)=vector_soft_row(C(:,idx)+Y2(:,idx)/u,wg(j)/u);
    end
    A=A-diag(diag(A));

    leq1=ones(1,n)*C-ones(1,n);
    leq2=C-A;
    stopC=max(max(abs(leq1(:))),max(abs(leq2(:))));
    
    
    if mod(i,1)==0 || stopC<tol
        disp(['itr ' num2str(i) ' ,stopALM=' num2str(stopC,'%2.3e')]);
%         plot(stopC(1:i)); pause(0.1);
    end
    if stopC<tol
        break;
    else
        Y1=Y1+u*leq1;
        Y2=Y2+u*leq2;
        u=min(u_max,ro*u);
    end
end