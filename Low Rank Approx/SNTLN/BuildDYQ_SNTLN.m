function DYQ = BuildDYQ_SNTLN(xk,m,n,k,alpha,theta)
% USED IN SNTLN function
%
% Construct Matrix DYQ, such that E_{k}(z)x = D^{-1}Y_{k}(x)Qz, where E_{k}(z) is a
% matrix of structured perturbations applied to S_{k}, where S_{k} = DTQ.
%
% Inputs
%
% xk : Vector
%
% m : Degree of polynomial f(x,y)
%
% n : Degree of polynomial g(x,y)
%
% k : Degree of GCD
%
% alpha : The optimal value of \alpha
%
% theta : The optimal value of \theta

% Get the number of coefficients in x1
nCoeffs_x1 = n-k+1;

% Split the vector x_{k} = [x_{1} x_{2}]
x1 = xk(1 : nCoeffs_x1);
x2 = xk(nCoeffs_x1+1:end);

% Build the matrix D^{-1}_{m+n-k}
D = BuildD_2Polys(m,n-k);

% Build the matrices Y_{1} and Y_{2}
Y1 = BuildT1(x1,m);
Y2 = BuildT1(x2,n);

% Build the matrices Q_{m} and Q_{n}
Qm = BuildQ1(m);
Qn = BuildQ1(n);

% Get a diagonal matrices of thetas corresponding to polynomials f(x) and
% g(x).
th_f = diag(theta.^(0:1:m));
th_g = diag(theta.^(0:1:n));

% Build DYQ
DYQ = D*[Y1*Qm*th_f  alpha.*Y2*Qn*th_g] ;


end






function Y1 = BuildDY_nchoosek_partition(xv,m,theta)
% Build the matrix D*T(xv)*R = D*T(f)*Q.
% Where R is a diagonal matrix of binomials and thetas corresponding to
% polynomial f(x). so f*R gives coefficients of f in the scaled modified
% Bernstein basis. a_{i}\theta^{i}\binom{m}{i}

global SETTINGS

switch SETTINGS.SYLVESTER_BUILD_METHOD
    case 'Standard'
        
        % Get the degree of xv
        n_t = GetDegree(xv);
        
        % Get the diagonal matrix D^{-1}
        D = BuildD_2Polys(n_t,m);
        
        % Get the binomials corresponding to f(x)
        Bi_m = GetBinomials(m);
        
        % Build the Cauchy matrix T(xv)
        T = BuildT1(xv,m);
        
        % Build the matrix R which is diagonal and includes thetas
        R = diag(GetWithThetas(Bi_m,theta));
        
        % Build the matrix DY
        Y1 = D*T*R;
        
    case 'Rearranged'
        
        %
        n_t = GetDegree(xv);
        
        %
        Y1 = zeros(n_t+1,m+1);
        
        
        for j = 0:1:m % for each column
            for i = j :1 : j + n_t 
                Y1(i+1,j+1) = ...
                    xv(i-j+1) .* theta^j ...
                    .* nchoosek(m+n_t-i,n_t-(i-j)) ...
                    .* nchoosek(i,i-j);
            end
        end
        
        %
        switch SETTINGS.BOOL_DENOM_SYL
            case 'y'
                Y1 = Y1./nchoosek(m+n_t,n_t);
            case 'n'
        end
end

end

function Y = BuildDY_log(xu,xv,t,alpha,theta)
% Construct Matrix Y, such that E_{k}x = Y_{k}z, where Ek is a matrix of
% structured perturbations.
%
%                           Inputs
%
%
% m - Degree of polynomial f
%
% n - Degree of polynomial g
%
% t - Degree of GCD
%
% mincol -
%
% x -
%
% theta -
%
%
%

% Global Variables
global SETTINGS

m_t = GetDegree(xu);
n_t = GetDegree(xv);

m = m_t + t;

n = n_t + t;

% First half
Y1 = zeros(m+n-t+1,m+1);

% Second half
Y2 = zeros(m+n-t+1,n+1);

% Build the first half of the sylvester matrix
% For each column k
for j=0:1:m
    % For each row i
    for i=j:1:j+length(xv)-1
        
        % Evaluate binomial coefficients in the numerator in logs
        Numerator_Eval_log = ...
            lnnchoosek(m+n-t-i,n-t-(i-j)) + ...
            lnnchoosek(i,j);
        
        % Convert to standard numeric form
        Numerator_Eval_exp = 10.^Numerator_Eval_log;
        
        
        Y1(i+1,j+1) = ...
            xv(i-j+1) .*(theta^j) .* Numerator_Eval_exp;
        
    end
end

% Build the second half of the sylvester matrix
for j=0:1:n
    for i=j:1:j+length(xu)-1
        
        % Evaluate the binomial coefficients in the numerator in logs
        Numerator_Eval_log = ...
            lnnchoosek(m+n-t-i,m-t-(i-j)) + ...
            lnnchoosek(i,j);
        
        % convert to standard numeric form
        Numerator_Eval_exp = 10.^Numerator_Eval_log;
        
        %
        Y2(i+1,j+1)= ...
            xu(i-j+1).*(theta^(j)) .* Numerator_Eval_exp;
    end
end



switch SETTINGS.BOOL_DENOM_SYL
    case 'y'
        % Include the denominator in the coefficient matrix.
        
        % Evaluate the binomial coefficient in the denominator of the first
        % partiton
        Denom1_eval_log = ...
            lnnchoosek(m+n-t,n-t);
        
        % Convert to standard numeric form.
        Denom1_eval_exp = 10.^Denom1_eval_log;
        
        % Divide the first partition by the constant common denominator.
        Y1 = Y1./Denom1_eval_exp;
        
        % Evaluate the binomial coefficient in the denominator of the second
        % partition.
        Denom2_eval_log = ...
            lnnchoosek(m+n-t,m-t);
        
        % Convert to standard numeric form
        Denom2_eval_exp = 10.^Denom2_eval_log;
        
        % Divide the second partition by the constant common denominator.
        Y2 = (alpha.*Y2)./Denom2_eval_exp;
        
    case 'n'
        % Exclude the denominator from the Sylvester Matrix
        Y2 = alpha.*Y2;
    otherwise
        error('SETTINGS.BOOL_DENOM_SYL is either y or n')
        
end

% Concatenate the partitions to form Y.
Y=[Y1,Y2];

end

