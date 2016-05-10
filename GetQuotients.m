function [ux,vx] = GetQuotients(fx_n,gx_n,t,alpha,theta)
% Given polynomials f(x) and g(x), get the quotient polynomials u(x) and
% v(x) such that f(x)*v(x) = g(x)*u(x).

global SETTINGS

% Get degree of input polynomial g(x)
n = GetDegree(gx_n);

% Get the polynomails f(\omega,\theta) and g(\omega,\theta)
fw = GetWithThetas(fx_n,theta);
gw = GetWithThetas(gx_n,theta);

% Build the t^th subresultant
St = BuildSubresultant(fw,alpha.*gw,t);

% Get the optimal column for removal
[opt_col] = GetOptimalColumn(St);

% Remove optimal column
At = St;
At(:,opt_col) = [];

% Get the optimal column c_{t} removed from S_{k}
ct = St(:,opt_col);

% Obtain the solution vector x = [-v;u]
x_ls = SolveAx_b(At,ct);

% Insert a zero into the position corresponding to the index of the optimal
% column so that S(f,g)*vec_x = 0.
vec_x =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ]  ;

% Obtain values for quotient polynomials u and v. still expressed in the
% scaled bernstein basis, including theta.
vw = vec_x(1:n-t+1);
uw = -vec_x(n-t+2:end);

% If Q is not included in Sylvester matrix, then binomials are included in
% x. Remove binomials

switch SETTINGS.BOOL_Q
    case 'y'
    case 'n'
        % Remove binomials from the coefficients.
        vw = GetWithoutBinomials(vw);
        uw = GetWithoutBinomials(uw);
        
    otherwise 
        error('err')
end

% Divide v(w) and u(w) to obtain u(x) and v(x)
vx = GetWithoutThetas(vw,theta);
ux = GetWithoutThetas(uw,theta);

end