function [S] = BuildSubresultant(fw,gw,k)
% BuildSubresultant(fw,gw,k,alpha)
% 
% This function builds the k-th subresultant matrix S_{k}, in the 
% Bernstein Basis.
%
%
%                       Inputs
%
% fx : Coefficients in Bernstein Basis
%
% gx : Coefficients in Bernstein Basis
%
% k:  Index of subresultant S_{k} to be built
%
% alpha :   Optimal value of alpha
%
%


%                       Global Variables.

global BOOL_Q
global SYLVESTER_BUILD_METHOD
%


% Get the degree of the polynomial f(x)
m = GetDegree(fw);

% Get the degree of the polynomial g.
n = GetDegree(gw);


switch SYLVESTER_BUILD_METHOD
    case 'Rearranged' % Build based on generating individual elements
        
        % Build First Partition   
        D1T1Q = BuildDT1Q1(fw,n-k);

        % Build Second Partition
        D2T2Q = BuildDT1Q1(gw,m-k);
        
        % Build Sylvester Matrix by concatenation of matrices A and B.
        S = [D1T1Q D2T2Q];
        
        
    case 'Standard' % Build based on multiplying by precalculated matrices .
        
        D = BuildD(m,n-k);
        T = BuildT(fw,gw,k);
        Q = BuildQ(m,n,k);

        switch BOOL_Q
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




