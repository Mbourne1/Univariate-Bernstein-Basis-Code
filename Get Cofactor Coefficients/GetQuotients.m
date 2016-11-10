function [ux,vx] = GetQuotients(fx,gx,t)
% Given polynomials f(x) and g(x), get the quotient polynomials u(x) and
% v(x) such that f(x)*v(x) = g(x)*u(x).

global SETTINGS

% Get degree of input polynomial g(x)
n = GetDegree(gx);

% Build the t^th subresultant
St = BuildSubresultant(fx,gx,t);

% Get the optimal column for removal
[~,idx_col] = GetMinDistance(St);

% Remove optimal column
At = St;
At(:,idx_col) = [];

% Get the optimal column c_{t} removed from S_{k}
ct = St(:,idx_col);

% Obtain the solution vector x = [-v;u]
x_ls = SolveAx_b(At,ct);

% Insert a zero into the position corresponding to the index of the optimal
% column so that S(f,g)*vec_x = 0.
vec_x =[
    x_ls(1:(idx_col)-1);
    -1;
    x_ls(idx_col:end);
    ]  ;

% Obtain values for quotient polynomials u and v. still expressed in the
% scaled bernstein basis, including theta.


switch SETTINGS.SYLVESTER_BUILD_METHOD
    case 'T'
        
        vx = vec_x(1:n-t+1);
        ux = -vec_x(n-t+2:end);
        
    case 'DT'
        
        vx_bi = vec_x(1:n-t+1);
        ux_bi = -vec_x(n-t+2:end);
        ux = GetWithoutBinomials(ux_bi);
        vx = GetWithoutBinomials(vx_bi);
        
    case 'DTQ'
        
        vx = vec_x(1:n-t+1);
        ux = -vec_x(n-t+2:end);
        
    case 'TQ'
        
        vx = vec_x(1:n-t+1);
        ux = -vec_x(n-t+2:end);
        
    case 'DTQ Rearranged Denom Removed'
        
        vx = vec_x(1:n-t+1);
        ux = -vec_x(n-t+2:end);
        
    case 'DTQ Rearranged'
        
        vx = vec_x(1:n-t+1);
        ux = -vec_x(n-t+2:end);
        
    otherwise 
        error('err')
end



end