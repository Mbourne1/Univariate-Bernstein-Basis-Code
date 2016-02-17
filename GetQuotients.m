function [ux,vx] = GetQuotients(fx_n,gx_n,t,alpha,theta)

% Get degrees of input polynomials
[r,c] = size(fx_n);
m = r-1;
[r,c] = size(gx_n);
n = r-1;


% Build the t^th subresultant
Sk = BuildSubresultant(fx_n,gx_n,t,alpha,theta);

% Get the optimal column for removal
[opt_col] = GetOptimalColumn(Sk);

% Remove optimal column
Aki = Sk;
Aki(:,opt_col) = [];
cki = Sk(:,opt_col);

inversionMethod = 0;

switch inversionMethod
    case 0 % Default use QR Decomposition
        [~,n2] = size(Aki);
        [Q,R] = qr(Aki);
        R1 = R(1:n2,:);
        cd = Q'*cki;
        c = cd(1:n2,:);
        x_ls = R1\c;
    case 1 % USE SVD for psuedoinverse
        x_ls = pinv(Aki)*cki;
end


% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ]  ;


% Obtain values for quotient polynomials u and v. still expressed in the
% scaled bernstein basis, including theta.
vw      =   vecx(1:n-t+1);
uw      =   -vecx(n-t+2:end);

% Divide vw and uw to obtain ux and vx

theta_nt = theta.^(0:1:n-t);
theta_mt = theta.^(0:1:m-t);

vx = vw./theta_nt';
ux = uw./theta_mt';

end