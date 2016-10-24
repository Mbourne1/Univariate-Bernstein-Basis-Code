function [S] = BuildSubresultant(fx,gx,k)
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


% Global Variables.
global SETTINGS


% Get the degree of the polynomial f(x)
m = GetDegree(fx);

% Get the degree of the polynomial g.
n = GetDegree(gx);


switch SETTINGS.SYLVESTER_BUILD_METHOD
    case 'Rearranged' % Build based on generating individual elements
        
        % Build First Partition   
        D1T1Q = BuildDT1Q1(fx,n-k);

        % Build Second Partition
        D2T2Q = BuildDT1Q1(gx,m-k);
        
        % Build Sylvester Matrix by concatenation of matrices A and B.
        S = [D1T1Q D2T2Q];
        
        
    case 'Standard' % Build based on multiplying by precalculated matrices .
        
        D = BuildD(m,n-k);
        T = BuildT(fx,gx,k);
        Q = BuildQ(m,n,k);

        switch SETTINGS.BOOL_Q
            case 'y' % Include Q
                S = D*T*Q;
            case 'n' % Exclude Q
                S = D*T;
            otherwise
                error('Error SETTINGS.BOOL_Q is either y or n')
        end
    otherwise
        error('SETTINGS.SYLVESTER_BUILD_METHOD is either standard or rearranged')
end

end




