function T = BuildT(fx,gx,alpha,theta,t)
%  Build the Toeplitz matrix T consisting of coefficients of fx and gx in
% the scaled bernstein basis. ie the entries of T are of the form
% a_{i}\binom{m}{i} and b_{i}\binom{n}{i}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                               Input

% fx: vector of coefficients of f(x) in the standard bernstein basis. a_{i}

% gx: vector of coefficients of g(x) in the standard Bernstein basis. b_{i}

% t : index of subresultant S_{t} to be formed. (Also degree of GCD)

%                              Output

% T : the partitioned matrix T = [T(f) T(g)].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get degree of polynomail f
m = length(fx) - 1;

% Get degree of polynomial g.
n = length(gx) - 1;

% Build Toeplitz matrix of f, the first partiton.
T1 = BuildT_Partition(fx,theta,n,t);

% Build Toeplitz matrix of g, the second partition.
T2 = BuildT_Partition(gx,theta,m,t);

% Concatenate the partitions.
T = [T1 alpha.*T2];
end

function T1 = BuildT_Partition(fx,theta,n,t)
% Build a Toeplitz Matrix of coefficients of fx.
% T1 \in \mathbb{R}^{(m+n-k+1)\times(n-k+1)}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs.


% fx : coefficients of polynomial f

% n : degree of polynomial g

% t :  index of subresultant S_{t} to be formed. (Also degree of GCD)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get degree of polynomail f
m = length(fx)-1;

% Initialise empty matrix T1, for storing Toeplitz T_{k}(f)
T1 = zeros(m+n-t+1,n-t+1);

% for each column of T1
for j = 0:1:n-t
    % for each coefficient f_{i} in fx
    for i=j:1:j+m
        T1(i+1,j+1) = fx(i-j+1) .*theta^(i-j) .* nchoosek(m,i-j);
    end
end



end