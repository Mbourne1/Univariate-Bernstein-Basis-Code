
function [F_max,F_min,G_max,G_min] = GetMaxMin(Sk,m,n,k)
% This function calculates the entries of minimum and maximum magnitude in
% the sylvester matrix S_{k} of each coefficient in the kth Sylvester
% subresulatnt matrix Sk.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Inputs.


% Sk  :  The modified Sylvester matrix of the polynomials F and G.
% m       :  The degree of F.

% n       :  The degree of G.

% k       :  The order of the subresultant matrix.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Outputs:


% F_max   :  A vector of length m+1, such that F_max(i) stores the
%            element of maximum magnitude of Sk that contains the
%            coefficient a(i) of F, i=1,...,m+1.

% F_min   :  A vector of length m+1, such that F_min(i) stores the
%            element of minimum magnitude of Sk that contains the
%            coefficient a(i) of F, i=1,...,m+1.

% G_max   :  A vector of length n+1, such that G_max(i) stores the
%            element of maximum magnitude of Sk that contains the
%            coefficient b(i) of G, i=1,...,n+1.

% G_min   :  A vector of length n+1, such that G_min(i) stores the
%            element of minimum magnitude of Sk that contains the
%            coefficient b(i) of G, i=1,...,n+1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sk=abs(Sk);

cols_f = n-k+1;   % the number of columns of F
cols_g = m-k+1;   % the number of columns of G

% Initialise the arrays that store the coefficients of minimum and
% maximum magnitude of F and G.
F_max = zeros(1,m+1);
F_min = zeros(1,m+1);
G_max = zeros(1,n+1);
G_min = zeros(1,n+1);

% Initialise this array for later.
F_sam = zeros(1,cols_f);
G_sam = zeros(1,cols_g);

% Consider initially the polynomial F. For each coefficient, select
% the largest and smallest elements that contain it.

% For each coefficient
for k1=1:1:m+1  
   % F_sam = zeros(1,cols_f);
    
    %   Put all the elements in the (k1-1)th leading diagonal, that is, the
    %   terms of f that contain the coefficient a(k1) in a vector and
    %   then locate the terms of minimum and maximum magnitude.
    for k2=1:1:cols_f
        F_sam(k2)=Sk(k1+k2-1,k2);
    end
    
    F_max(k1)=max(abs(F_sam));
    F_min(k1)=min(abs(F_sam));
    
end

% Repeat this procedure for the coefficients b(i) of G.


for k1=1:1:n+1
    
    for k2=cols_f+1:1:cols_f+cols_g
        G_sam(k2-cols_f)=Sk(k1+k2-1-cols_f,k2);
    end
    
    G_max(k1)=max(abs(G_sam));
    G_min(k1)=min(abs(G_sam));
    
end

