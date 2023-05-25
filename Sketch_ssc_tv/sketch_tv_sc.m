function Z=sketch_tv_sc(X,B,opts)
% X: input matrix, each column of X is a data point
% B: sketched dictionary with X*R, R is a random matrix
% opts: related parameters

lambda=opts.lam;
lambda_tv=opts.lam_tv;
tol=opts.tol;
rows=opts.row;
cols=opts.col;
Max_itr=opts.max_itr;
ro=opts.ro;
step=opts.step;

n=size(B,2);
N=rows*cols;
Z=zeros(n,N);
A=zeros(n,N);

u1=zeros(rows,cols,n);
u2=zeros(rows,cols,n);
Y1=zeros(n,N);
Y2=zeros(n,N);
Y31=zeros(rows,cols,n);
Y32=zeros(rows,cols,n);
Btx=B'*X;
Btb=B'*B;


% define total variation operators (compute by fft and ifft)
eigDtD      = abs(fftn([1 -1],  [rows cols n])).^2 + abs(fftn([1 -1]', [rows cols n])).^2;
[D,Dt]      = defDDt();

eigA  = 2 + eigDtD;


for itr=1:Max_itr
    % solve f-subproblem
    C=(Btb+ro*eye(n))\(Btx+ro*A+Y1);
    tp=Z+C-(Y1+Y2)/ro;
    rhs=fftn(reshape(tp',rows,cols,n)+Dt(u1-Y31/ro,  u2-Y32/ro));
    A_cub     = real(ifftn(rhs./eigA));
    A=reshape(A_cub,rows*cols,n)';
    Z=soft(A+Y2/ro,lambda/ro);
    
    % solve u-subproblem
    [Df1,Df2] = D(A_cub);
    v1 = Df1+Y31/ro;
    v2 = Df2+Y32/ro;
    u1= soft(v1,lambda_tv/ro);
    u2= soft(v2,lambda_tv/ro);
    
    AC=A-C;
    AZ=A-Z;
    Dd1=Df1-u1;
    Dd2=Df2-u2;
    Y1=Y1+ro*AC;
    Y2=Y2+ro*AZ;
    Y31=Y31+ro*Dd1;
    Y32=Y32+ro*Dd2;
    
    err(itr)=max([max(AC(:)) max(AZ(:)) max(Dd1(:)) max(Dd2(:))]);
    if mod(itr,20)==0
        fprintf('%3g \t %6.4e \n',itr,err(itr));
    end
    if err(itr)<tol
        break
    else
        ro=ro*step;
    end
end

 end

function [D,Dt] = defDDt()
D  = @(U) ForwardD(U);
Dt = @(X,Y) Dive(X,Y);
end

function [Dux,Duy] = ForwardD(U)
Dux = [diff(U,1,2), U(:,1,:) - U(:,end,:)];% difference in data cube, output is a data cube
Duy = [diff(U,1,1); U(1,:,:) - U(end,:,:)];
end

function DtXY = Dive(X,Y)
DtXY = [X(:,end,:) - X(:, 1,:), -diff(X,1,2)];
DtXY = DtXY + [Y(end,:,:) - Y(1, :,:); -diff(Y,1,1)];
end