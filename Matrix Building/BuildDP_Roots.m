
function [P] = BuildDP_Roots(m,n,alpha,theta,q,t,ratio)
% See [Report - SNTLN - Roots - Derivation of the first rearrangement]
%
%
%                       Inputs
%
% m :   Degree of polynomial f
%
% n :   Degree of polynomial g
%
% alpha :   Optimal value of alpha
%
% theta :   Optimal value of theta
%
% idxMinCol :   Index of optimal column
%
% t :   Degree of GCD and index of subresultant S_{t}
%
% ratio :   Ratio of geometric means \frac{\lamdba}{\mu}
%
%

% where q is the index of the column removed from the Sylvester matrix. q =
% 1,...,m+n-2t+2

if q <= n-t+1

    P = buildP_LHS_Roots(m,n,theta,q,t);
    
elseif q <= m+n-(2*t)+2

    P = buildP_RHS_Roots(m,n,alpha,theta,q,t,ratio);
end




end


function [P] = buildP_LHS_Roots(m,n,theta,idxMinCol,t)
% See [Report - SNTLN - Roots - Derivation of the first rearrangement]
%
%
%
% Inputs
%
% m :   Degree of polynomial f(x)
%
% n :   Degree of polynomial g(x)
%
% theta :   Optimal value of theta
%
% idxMinCol :   Index of optimal column
%
% t :   Degree of GCD d(x) and index of subresultant S_{t}


qhat = idxMinCol-1;

Z1 = zeros(qhat,m+1);

x = (m+n-t+1) - qhat - (m+1);
Z2 = zeros(x,m+1);

G = zeros(m+1,m+1);
for i = 0:1:m
    G(i+1,i+1) = 1 .* (theta^i) .* nchoosek(i+qhat,qhat) .* nchoosek(m+n-t-(i+qhat),m-i) ./ nchoosek(m+n-t,n-t);
end


P = [
    Z1;
    G; 
    Z2];



end

function [P] = buildP_RHS_Roots(m,n,alpha, theta, idxMinCol,t,ratio)
% See [Report - SNTLN - Roots - Derivation of the first rearrangement]
%
%
%                       Inputs
%
% m :   Degree of polynomial f
%
% n :   Degree of polynomial g
%
% alpha :   Optimal value of alpha
%
% theta :   Optimal value of theta
%
% idxMinCol :   Index of optimal column
%
% t :   Degree of GCD and index of subresultant S_{t}
%
% ratio :   Ratio of geometric means \frac{\lamdba}{\mu}

% let j be the index of the column which has been removed from the right
% hand partition, from 1,...,m-t

j = idxMinCol - (n-t+1);

jhat = j-1;

Z1 = zeros(jhat,m+1);

y = (m+n-t+1)- jhat - m;

Z2 = zeros(y,m+1);

G = zeros(m,m+1);

% for each of the coefficients of c_{q}
for i = 0:1:m-1

    G(i+1,i+1) = -  m .* (theta^i) .* nchoosek(i+jhat,jhat) .* nchoosek(m+n-t-(i+jhat),n-i) ./ nchoosek(m+n-t,m-t);
    G(i+1,i+2) =    m .* (theta^i) .* nchoosek(i+jhat,jhat) .* nchoosek(m+n-t-(i+jhat),n-i) ./ nchoosek(m+n-t,m-t);
    
    
end


P = [Z1;
    alpha.*ratio.*G;
    Z2];
end


