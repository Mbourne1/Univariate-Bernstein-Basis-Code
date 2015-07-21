function C1 = BuildC1(uw,t)
% Used in APF, build the matrix C1_{t}(u) or C1_{t}(v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs


%   uw :    Polynomial u including 

%   t : Degree of the gcd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Outputs


% C1 :  The Toeplitz structured matrix containing coefficients of
%       polynomial u in with its binomial coefficients included.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get degree of polynomial u
m_t = length(uw)-1;

C1 = zeros(m_t+t,t);



% for each column in C1 0,...,t

for j = 0:1:t
    % for each coefficient u_{0}...,u_{m-k}
    for i = j:1:m_t+j
        C1(i+1,j+1) = uw(i-j+1) .* nchoosek(m_t,i-j);
    end
end

end