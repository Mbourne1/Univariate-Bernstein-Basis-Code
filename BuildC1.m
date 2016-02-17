function C1 = BuildC1(uw,t)
% Used in APF, build the matrix C1_{t}(u) or C1_{t}(v)
%
%                           Inputs
%
%
%   uw :    Polynomial u(w) 
%
%   t :     Degree of the gcd
%
%
%                           Outputs
%
%
% C1 :  The Toeplitz structured matrix containing coefficients of
%       polynomial u in with its binomial coefficients included.
%
%

% Get degree of polynomial u
[r,c] = size(uw);
m_t = r -1;

% Initialise the matrix
C1 = zeros(m_t+t,t);

% Include binomial coefficients $u(w)_{i}.*nchoosek(m-t,i)$
uw_bi = uw .* GetBinomials(m_t);

% for each column in C1 0,...,t
for j = 0:1:t
    % Insert the coefficients of u(w) into the column.
    C1(j+1:m_t+1+j,j+1) = uw_bi;
end

end