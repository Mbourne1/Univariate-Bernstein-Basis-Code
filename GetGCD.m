function [dw] = GetGCD(uw,vw,fw_n,gw_n,t)
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
m = length(fw_n)-1;

% Get degree of polynomial g
n = length(gw_n)-1;

% Build solution vector bk = [f;g]
bk = [fw_n ; gw_n];

% Build the coefficient vector HCG

switch Bool_APFBuildMethod
    case 1 % use my method.
        switch bool_log
            case 1 % use logs
                C1 =  BuildC1_log(uw,m,t);
                C2 =  BuildC1_log(vw,n,t);
                
            case 0 % use nchoosek
                C1 =  BuildC1_nchoosek(uw,m,t);
                C2 =  BuildC1_nchoosek(vw,n,t);
        end
        HCG = [C1 ; C2];
    case 0 % use premade matrix method
        
        H = BuildH(m,n);
        C = BuildC(uw,vw,t);
        G = BuildG(t);
        
        HCG = H*C*G;
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
    case 1 % Denominator included
                
        DenomEval_log = lnnchoosek(m,t);
        
        DenomEval_exp = 10.^DenomEval_log;
        
        C = C./DenomEval_exp;
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
    case 1
        C = C./nchoosek(m,t);
end

end


function H = BuildH(m,n)
% Build the diagonal matrix H_{k}^{-1} consisting of binomial coefficients
% corresponding to f and g. The Coefficient matrix is given by
% H_{k}^{-1}C_{k}(u,v)G_{k}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% m : degree of polynomial f
% n : degree of polynomial g
% t : degree of gcd.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build partition H1, diagonal binomial coefficients corresponding to
% polynomial f.
H1 = BuildH_Partition(m);

% Build partition H2, diagonal binomial coefficients corresponding to
% polynomial g.
H2 = BuildH_Partition(n);

% Form diagonal matrix of the two partitions.
H = blkdiag(H1,H2);
end

function H1 = BuildH_Partition(m)
%% Build partition of the Matrix H.

% Initalise empty vector H1
H1 = zeros(1,m+1);

% For each element in H1, get the binomial coefficient \binom{m}{i}
for i=0:1:m
    H1(i+1) = 1./nchoosek(m,i);
end

% Form matrix from vector by diagonalizing H1.
H1 = diag(H1);

end

function C = BuildC(u,v,t)
%% Build the Matric C(f,g) = [C(f) | C(g)] \in \mathbb{R}^{(m+n+2)\times(k+1)}

% Build partition of C corresponding to polynomial u
C1 = BuildC_Partition(u,t);

% Build partition of C corresponding to polynomial v
C2 = BuildC_Partition(v,t);

C = [C1 ; C2];
end

function C1 = BuildC_Partition(u,t)

m = length(u)-1 +t;

C1 = zeros(m+1,t+1);

% for each column in C1
for j=0:1:t
    % for each coefficient u_{i} in u
    for i=j:1:j+m-t
        C1(i+1,j+1) = u(i-j+1) .* nchoosek(m-t,i-j);
    end
end
end

function G = BuildG(t)

G = zeros(1,t+1);

for i=0:1:t
    G(i+1) = nchoosek(t,i);
end


G = diag(G);
end
