function [ux,vx] = GetQuotients(fx_n,gx_n,t,alpha,theta)
% Given polynomials f(x) and g(x), get the quotient polynomials u(x) and
% v(x) such that f(x)*v(x) = g(x)*u(x).

global BOOL_Q

% Get degrees of input polynomial f(x)
[r,~] = size(fx_n);
m = r-1;

% Get degree of input polynomial g(x)
[r,~] = size(gx_n);
n = r-1;

% Get the polynomails f(\omega,\theta) and g(\omega,\theta)
fw = fx_n .* (theta.^(0:1:m)');
gw = gx_n .* (theta.^(0:1:n)');

% Build the t^th subresultant
St = BuildSubresultant(fw,gw,t,alpha);

% Get the optimal column for removal
[opt_col] = GetOptimalColumn(St);

% Remove optimal column
At = St;
At(:,opt_col) = [];
ct = St(:,opt_col);

x_ls = SolveAx_b(At,ct)


% Obtain the solution vector x = [-v;u]
vec_x =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ]  ;


% Obtain values for quotient polynomials u and v. still expressed in the
% scaled bernstein basis, including theta.
vw      =   vec_x(1:n-t+1);
uw      =   -vec_x(n-t+2:end);

% If Q is not included in Sylvester matrix, then binomials are included in
% x. Remove binomials

switch BOOL_Q
    case 'y'
    case 'n'
        Bi_mt = GetBinomials(m-t);
        Bi_nt = GetBinomials(n-t);
        
        vw = vw./ Bi_nt;
        uw = uw./ Bi_mt;
    otherwise 
        error('err')
end

% Divide v(w) and u(w) to obtain u(x) and v(x)
theta_nt = theta.^(0:1:n-t);
theta_mt = theta.^(0:1:m-t);

vx = vw./theta_nt';
ux = uw./theta_mt';

end