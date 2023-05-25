% function NMI = getNMI(A,B)
% % call: 
% % 
% %         NMI = getNMI(A,B)
% %         
% % This function computes the Normalized Mutual Information (NMI) between 
% % 2 modular partitions or community structures given in vectors A and B.
% % NMI is a measure of the similarity between the two graph partitions, and 
% % its interpretation follows that of canonical Mutual Information, i.e. 
% % assuming that I have complete knowledge of partition A, how much knowledge 
% % can I get, from A, about partition B? If A and B are identical, then 
% % NMI = 1, whereas if A and B are independent, NMI --> 0.
% % 
% % INPUT
% % 
% %      A    :   Community partition of graph A *
% %      B    :   Community partition of graph B *
% % 
% %      * A and B are N-length vectors, where each i-th element is the integer 
% %        labeling the k-th community to which node i-th was assigned.
% %        
% %        
% % OUTPUT
% % 
% %      NMI  :   Normalized Mutual Information 
% %      
% %      
% %     
% % References:
% % [1] Kuncheva & Hadjitodorov, 2004, IEEE Intern. Conf. on Syst. Man and Cybern.
% % [2] Alexander-Bloch et al., 2012, NeuroImage
% % 
% % 
% % R. G. Bettinardi
% % ------------------------------------------------------------------------------------------------------------------------------------------------------------
% if nargin < 2
%     error('One of the two inputs is missing !!!')
% end
% % ASSURE A; B INTEGERS!        
% N  = numel(A);
% Ca = max(A);
% Cb = max(B);
% N1 = zeros(Ca*Cb,1);
% D1 = zeros(Ca,1);
% D2 = zeros(Cb,1);
% n1_k = 0;
% for i = 1:Ca
%     
%     N_idot   = sum(A==i);
%     D1(i) = N_idot * log(N_idot/N);
%     
%    for j = 1:Cb 
%        
%        n1_k = n1_k + 1;
%        
%        N_dotj = sum(B==j);
%        N_ij   = sum(((A==i) + (B==j))==2);
%        
%        if N_ij == 0
%            N1(n1_k) = 0;
%        else
%            N1(n1_k) = N_ij * log( (N_ij*N)/(N_idot*N_dotj)  );
%        end
%        
%        D2(j) = N_dotj * log(N_dotj/N);
%        
%    end
% end
% NMI = (-2 * sum(N1) ) / ( sum(D1) + sum(D2) );

%% another implementation
function z = getNMI(x, y)
% Compute normalized mutual information I(x,y)/sqrt(H(x)*H(y)) of two discrete variables x and y.
% Input:
%   x, y: two integer vector of the same length 
% Ouput:
%   z: normalized mutual information z=I(x,y)/sqrt(H(x)*H(y))
% Written by Mo Chen (sth4nth@gmail.com).
assert(numel(x) == numel(y));
n = numel(x);
x = reshape(x,1,n);
y = reshape(y,1,n);

l = min(min(x),min(y));
x = x-l+1;
y = y-l+1;
k = max(max(x),max(y));

idx = 1:n;
Mx = sparse(idx,x,1,n,k,n);
My = sparse(idx,y,1,n,k,n);
Pxy = nonzeros(Mx'*My/n); %joint distribution of x and y
Hxy = -dot(Pxy,log2(Pxy));


% hacking, to elimative the 0log0 issue
Px = nonzeros(mean(Mx,1));
Py = nonzeros(mean(My,1));

% entropy of Py and Px
Hx = -dot(Px,log2(Px));
Hy = -dot(Py,log2(Py));

% mutual information
MI = Hx + Hy - Hxy;

% normalized mutual information
z = sqrt((MI/Hx)*(MI/Hy));
z = max(0,z);