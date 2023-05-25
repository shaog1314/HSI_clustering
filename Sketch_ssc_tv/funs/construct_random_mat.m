function [R] = construct_random_mat(N,M,random_type,silent)
%CONSTRUCT_RANDOM_MAT Function that generates random JLT matrices.
%   Inputs: N,M - dimensions of the random matrix
%   random_type - The type of JLT matrix to be constructed, i.e. 'Normal',
%   'Rademacher', 'Sparse', etc.
%        silent - whether or not to show prompts
%   Output:   R - NxM JLT matrix of type random_type 
%   Panagiotis Traganitis traga003@umn.edu

if nargin < 3
    random_type = 'Normal';
end
if nargin < 4
    silent = false;
end
R = zeros(N,M);
if M < N
    const = 1/sqrt(M);
else
    const = 1/sqrt(N);
end

switch lower(random_type)
    case {'normal'}
        if ~silent
            disp('Generating random matrix with i.i.d. Normal entries');
        end
        R = const.*randn(N,M);
    case {'rademacher'}
        if ~silent
        disp('Generating random matrix with i.i.d. Rademacher entries');
        end
        R = randi([0 1],N,M);
        R(~R(:)) = -1; R = const.*R;
    case {'sparse'}
        if ~silent
            disp('Generating random sparse embedding matrix');
        end
        indx = randi([1 M],N,1);
        Phi = sparse(indx,(1:N)',1,M,N);
        val = randi([0 1],N,1); val(val == 0) = -1;
        D = spdiags(val,0,N,N);
        R = Phi*D; R = R';
    case {'countsketch'}
        if ~silent
        disp('Generating random CountSketch matrix');
        end
        indx = randi([1 N],M,1); val = randi([0 1],M,1); val(val == 0) = -1;
        Rtmp = [indx (1:M)' val];
        R = spconvert(Rtmp);
        if size(R,1) < N
            R = [R; zeros(N - size(R,1),M)];
        end     
    case {'hadamard'}
        if ~silent
        disp('Generating ROS Hadamard matrix');
        end
        tt = log2(N);
        if rem(tt,2) == 0 
            H = hadamard(N); D = randi([0 1],N,1); D(D == 0) = -1; D = diag(D); I = eye(N); indx = randperm(N);
            R = sqrt(N/M)*D*H'*I(:,indx(1:M));
        else
            T = ceil(tt); New_N = 2^T;
            H = hadamard(New_N);  D = randi([0 1],New_N,1); D(D == 0) = -1; D = diag(D); I = eye(New_N); indx = randperm(New_N);
            R = sqrt(New_N/M)*D*H'*I(:,indx(1:M));
            %error('N is not a power of 2 - cannot create Hadamard matrix');
        end
    otherwise
        disp('No valid syntax provided - returning all zeros matrix');
end
end

