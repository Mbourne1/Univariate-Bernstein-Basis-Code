function DY = BuildDY(m,n,t,opt_col,x,alpha,theta)
% USED IN SNTLN function

% Construct Matrix DY, such that E_{k}(z)x = D^{-1}Y_{k}(x)z, where E_{k}(z) is a 
% matrix of structured perturbations applied to S_{k}, where S_{k} = DTQ.
% Similarly Y_{k} can be used in the following
%   S_{k}x = Y_{k}[f;g]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs

% m - Degree of polynomial f

% n - Degree of polynomial g

% t - Degree of GCD

% opt_col - the optimal column (also referred to as ck) for removal from 
% the Sylvester matrix such that the residual obtained is minimal. 
% Ak x = ck.

% x - the Approximate solution, which consists of coefficients of the
% coprime polynomials u and v.

% alpha - the calculated optimal value of alpha

% theta - the calculated optimal value of theta


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                            Global Variables


% bool_log - (Boolean)
%   1 :- Perform calculations by log method
%   0 :- Perform calculations by standard method.

global bool_log

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch bool_log
    case 'n' % Use nchoosek method
        
        DY = BuildDY_nchoosek(m,n,t,opt_col,x,alpha,theta);
        
    case 'y' % Use Logs method
        
        DY = BuildDY_log(m,n,t,opt_col,x,alpha,theta);
end
end


function DY = BuildDY_nchoosek(m,n,t,mincol,x_ls_wrt_w,alpha,theta)

% Construct Matrix Y, such that E_{k}x = Y_{k}z, where Ek is a matrix of
% structured perturbations applied to S_{k} = DTQ

% Build in the rearranged format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                               Inputs


% m - Degree of polynomial f

% n - Degree of polynomial g

% t - Degree of GCD

% mincol -

% theta -

% alpha -


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Global Variables

% bool_denom_syl - (Boolean) Given the rearrangement of the Sylvester matrix in
% the Bernstein basis, each partition of each subresultant has a common
% divisor to its elements.
%    1 :- Include Common Denominator.
%    0 :- Exclude Common Denominator.
global bool_denom_syl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xa = x_ls_wrt_w(1:mincol-1) ;
xb = x_ls_wrt_w(mincol:end) ;
x_ls_wrt_w = [xa; 0 ;xb] ;% Insert zero into vector

% Get number of cols in left partition
c_leftmatrix = (n-t+1);

% Get number of cols in right partition
c_rightmatrix = (m-t+1);


Xv_w = x_ls_wrt_w(1:c_leftmatrix);

Xu_w = x_ls_wrt_w(c_leftmatrix+1:c_leftmatrix+c_rightmatrix);


%Build Empty Sylvester Matrix
% First half
Y1 = zeros(m+n-t+1,m+1);

% Second half
Y2 = zeros(m+n-t+1,n+1);

% Build the first half of the matrix DY

% For each column j = 0,...,
for j=0:1:m
    % For each row i = j,...,
    for i=j:1:j + length(Xv_w)-1
        Y1(i+1,j+1) = ...
            Xv_w(i-j+1) .*(theta^(j)) ...
            .* nchoosek(m+n-t-i,n-t-(i-j)) ...
            .* nchoosek(i,i-j);
    end
end


% Build the second half of the matrix DY
% for each column j = 0,...,n
for j=0:1:n
    % for each row i = j,...,
    for i=j:1:j + length(Xu_w)-1
        Y2(i+1,j+1)= ...
            Xu_w(i-j+1).*(theta^(j)) .* ...
            nchoosek(m+n-t-i,m-t-(i-j)) .*...
            nchoosek(i,i-j);
    end
%X2(i-j+1).*(theta^(j)) .* ...
end


switch bool_denom_syl
    case 'y'
        %Include the denominator
        Y1 = Y1./nchoosek(m+n-t,n-t);
        Y2 = alpha.*Y2./nchoosek(m+n-t,m-t);
    case 'n'
        % Exclude the denominator
        
        Y2 = alpha.*Y2;
end


DY=[Y1,Y2];

end

function Y = BuildDY_log(m,n,t,mincol,x,alpha,theta)
% Construct Matrix Y, such that E_{k}x = Y_{k}z, where Ek is a matrix of
% structured perturbations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs


% m - Degree of polynomial f

% n - Degree of polynomial g

% t - Degree of GCD

% mincol -

% x -

% theta -


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                        Global Variables

% BOOL_DENOM - (Boolean) Given the rearrangement of the Sylvester matrix in
% the Bernstein basis, each partition of each subresultant has a common
% divisor to its elements.
%    1 :- Include Common Denominator.
%    0 :- Exclude Common Denominator.
global bool_denom_syl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xa = x(1:mincol-1) ;
xb = x(mincol:end) ;
x = [xa; 0 ;xb] ;% Insert zero into vector

c_leftmatrix = (n-t+1);
c_rightmatrix = (m-t+1);

X1 = x(1:c_leftmatrix);
X2 = x(c_leftmatrix+1:c_leftmatrix+c_rightmatrix);

%Build Empty Sylvester Matrix

% First half
Y1 = zeros(m+n-t+1,m+1);

% Second half
Y2 = zeros(m+n-t+1,n+1);

% Build the first half of the sylvester matrix
% For each column k
for j=0:1:m
    % For each row i
    for i=j:1:j+length(X1)-1
        
        % Evaluate binomial coefficients in the numerator in logs
        Numerator_Eval_log = ...
            lnnchoosek(m+n-t-i,n-t-(i-j)) + ...
            lnnchoosek(i,j);
        
        % Convert to standard numeric form
        Numerator_Eval_exp = 10.^Numerator_Eval_log;
        
        
        Y1(i+1,j+1) = ...
            X1(i-j+1) .*(theta^j) .* Numerator_Eval_exp;
        
    end
end

% Build the second half of the sylvester matrix
for j=0:1:n
    for i=j:1:j+length(X2)-1
        
        % Evaluate the binomial coefficients in the numerator in logs
        Numerator_Eval_log = ...
            lnnchoosek(m+n-t-i,m-t-(i-j)) + ...
            lnnchoosek(i,j);
        
        % convert to standard numeric form
        Numerator_Eval_exp = 10.^Numerator_Eval_log;
        
        %
        Y2(i+1,j+1)= ...
            X2(i-j+1).*(theta^(j)) .* Numerator_Eval_exp;
    end
end



switch bool_denom_syl
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
        
end

% Concatenate the partitions to form Y.
Y=[Y1,Y2];

end

