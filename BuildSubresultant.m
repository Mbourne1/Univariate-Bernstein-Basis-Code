function [S] = BuildSubresultant(fx,gx,k,alpha,theta)
% This function calculates the kth subresultant matrix S_{k}, in the 
% Bernstein Basis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Inputs
%
% fx: Coefficients in Bernstein Basis
%
% gx: Coefficients in Bernstein Basis
%
% theta :   Optimal value of theta
%
% alpha :   Optimal value of alpha
%
%   k:  Index of subresultant S_{k} to be built
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Variables.

global bool_q
global bool_sylvesterBuildMethod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get the degree of the polynomial f. 
[r,~] = size(fx);
m = r - 1;

% Get the degree of the polynomial g.
[r,~] = size(gx);
n = r - 1;



switch bool_sylvesterBuildMethod
    case 'rearranged' % Build based on generating individual elements
        
        % Build First Partition   
        A = BuildDT1Q1(fx,theta,n,k);
        
        % Build Second Partition
        B = BuildDT1Q1(gx,theta,m,k);
        
        % Build Sylvester Matrix by concatenation of matrices A and B.
        S = [A alpha.*B];
        
        
    case 'standard' % Build based on multiplying by precalculated matrices .
        
        D = BuildD(m,n,k);
        T = BuildT(fx,gx,alpha,theta,k);
        Q = BuildQ(m,n,k);

        switch bool_q
            case 'y' % Include Q
                S = D*T*Q;
            case 'n' % Exclude Q
                S = D*T;
            otherwise
                error('err')
        end
    otherwise
        error('bool_sylvesterBuildMethod is either standard or rearranged')
end

end




