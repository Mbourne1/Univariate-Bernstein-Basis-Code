function [dx] = GetGCD(ux,vx,fx_n,gx_n,t,alpha,theta)
% Get The Coefficients of the approximate GCD using Quotient Polynomials.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Inputs

% uw = quotient of f where uw is in the form u_{i}\theta^{i}.

% vw = quotient of g where vw is in the form v_{i}\theta^{i}.

% fw_n = coefficients of polynomial f in modified bernstein basis.

% gw_n = coefficients of polynomial g in modified bernstein basis.

% t = degree of AGCD.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Global Variables

global bool_log

global Bool_APFBuildMethod


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get degree of polynomials f
[r,c] = size(fx_n);
m = r-1;
% Get degree of polynomial g
[r,c] = size(gx_n);
n = r-1;

% Get fw and gw
fw_n = fx_n .* theta.^(0:1:m)';
gw_n = gx_n .* theta.^(0:1:n)';

% Get uw and vw
uw = ux .* theta.^(0:1:m-t)';
vw = vx .* theta.^(0:1:n-t)';

% Build solution vector bk = [f;g]
bk = [fw_n ; alpha .* gw_n];

% Build the coefficient vector HCG

switch Bool_APFBuildMethod
    case 'rearranged' % use rearranged method with common denominator
        switch bool_log
            case 'y' % use logs
                C1 =  BuildC1_log(uw,m,t);
                C2 =  BuildC1_log(vw,n,t);
                
            case 'n' % use nchoosek
                C1 =  BuildC1_nchoosek(uw,m,t);
                C2 =  BuildC1_nchoosek(vw,n,t);
        end
        HCG = [C1 ; C2];
    case 'standard' % use premade matrix method
        
        H = BuildH(m,n);
        C = BuildC(uw,vw,t);
        G = BuildG(t);
        
        HCG = H*C*G;
    otherwise
        error('Bool_APFBuildMethod must either be (standard) or (rearranged)')
end


try
    
    inversionMethod = 1;
    switch inversionMethod
        case 0 % Default use QR decomposition
            
            [~,n2] = size(HCG);
            [Q,R] = qr(HCG);
            R1 = R(1:n2,:);
            cd = Q'*bk;
            c = cd(1:n2,:);
            x_ls = R1\c;
            
            dw = x_ls;
        case 1 % Use SVD for psuedoinverse
            
            x_ls = pinv(HCG)*bk;
            dw = x_ls;
            
            
    end
    
catch
    fprintf('Could not perform QR Decomposition used in obtaining the coefficients of the GCD.')
end

dx = dw./(theta.^(0:1:t)');


end

function [C] = BuildC1_log(uw,m,t)
%% Build The Toeplitz Matrix of Quotient Polynomial C_{1}(u)
%% where C_{1} \in \mathbb{R}^{(m+1)\times(t+1)}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   u : coefficients of quotient polynomial in bernstein basis
%   t : degree of gcd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global bool_denom_apf

% Initialise Toeplitz Matrix.
C = zeros(m+1,t+1);

% for each column of the matrix
for j= 0:1:t
    % for each element u_{i} of u(w)
    for i = j:1:j+(m-t)
        
        % Evaluate the binomial coefficient in the numerator by log method.
        NumEval_log = lnnchoosek(m-i,t-j) + lnnchoosek(i,j);
        
        % Get exponent of the log of the numrator binomial coefficient.
        NumEval_exp = 10.^NumEval_log;
        
        % Multiply coefficient u_{i-j} by binomial coefficient.
        C(i+1,j+1) = uw(i-j+1) .* NumEval_exp ;
        
    end
    
end


switch bool_denom_apf
    case 'y' 
        % Denominator included
                
        DenomEval_log = lnnchoosek(m,t);
        
        DenomEval_exp = 10.^DenomEval_log;
        
        C = C./DenomEval_exp;
    case 'n'
        % Exclude Denominator
    otherwise 
        error('bool_denom_apf')
end

end


function [C] = BuildC1_nchoosek(uw,m,t)
global bool_denom_apf

C = zeros(m+1,t+1);

for j= 0:1:t
    for i = j:1:j+(m-t)
        C(i+1,j+1) = uw(i-j+1) .* nchoosek(m-i,t-j) .* nchoosek(i,j);
    end
end

switch bool_denom_apf
    case 'y'
        % Include the denominator
        C = C./nchoosek(m,t);
    case 'n'
        % Exclude the denominator
    otherwise
        error('bool_denom_apf must be either y or n')
end

end




