function [S] = Subresultant_BernsteinBasis(fx,gx,theta,alpha,k)
% This function calculates the kth subresultant matrix S_{k}, in the Bernstein Basis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Inputs

% fx: Coefficients in Bernstein Basis

% gx: Coefficients in Bernstein Basis

% theta :   Optimal value of theta

% alpha :   Optimal value of alpha

%   k:  Index of subresultant S_{k} to be built


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Variables.

% bool_q (Boolean)
%   1 :- Q included in the Sylvester Matrix S(f,g) = D^{-1}T(f,g)Q
%   0 :- Q excluded from Sylvester Matrix S(f,g) = D^{-1}T(f,g)
global bool_q

% Bool_sylvesterBuildMethod
% 1 :   Build based on individual elements of the Sylvester matrix, each
%       (i,j) element is calculated independently.
% 0 :   Use Naive method calculate D, calculate S, calculate Q, then
%       calculate DTQ.
global bool_sylvesterBuildMethod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get the degree of the polynomial f. 
% Get the degree of the polynomial g.
m = length(fx) - 1;
n = length(gx) - 1;



switch bool_sylvesterBuildMethod
    case 1 % Build based on generating individual elements
        
        % Build first partition
        A = BuildToeplitz(fx,theta,n,k);
        
        % Build Second Partition
        B = BuildToeplitz(gx,theta,m,k);
        
        % Build Sylvester Matrix by concatenation of matrices A and B.
        S = [A alpha.*B];
        
        
    case 0 % Build based on multiplying by precalculated matrices .
        
        D = BuildD(m,n,k);
        T = BuildT(fx,gx,alpha,theta,k);
        Q = BuildQ(m,n,k);

        switch bool_q
            case 1 % Include Q
                S = D*T*Q;
            case 0 % Exclude Q
                S = D*T;
        end

end

end




