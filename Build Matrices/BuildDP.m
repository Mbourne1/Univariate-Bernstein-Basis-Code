function [P] = BuildDP(m,n,theta,mincol,t)
% BuildDP(m,n,theta,mincol,t)
%
% Build the matrix DP.
%
% Used in SNTLN.m
%
% Build the matrix DP such that DP * z gives the column of the Sylvester
% matrix whose column index is equal to mincol
%
%
% Inputs
%
% m :   Degree of polynomial f(x)
%
% n :   Degree of polynomial g(x)
%
% theta :   Optimal value of theta
%
% num_mincol :  Index of column c_{k} removed from S_{k}(f,g), index starts
%               at 1
%
% t :   Degree of GCD d(x)
%
% Outputs.
%
% P : matrix DP


if mincol <= (n-t+1) % the optimal column is from poly f
    
    r_aboveG = mincol-1;
    r_belowG = n-t+1-mincol;
    
    p1 = zeros(r_aboveG,m+1);
    p2 = zeros(r_aboveG,n+1);
    
    G  = LeftDG_Build(m,n,t,mincol,theta);
    p3 = zeros(m+1,n+1);
    
    p5 = zeros(r_belowG,m+1);
    p6 = zeros(r_belowG,n+1);
    
    P=[ p1  p2;...
        G   p3;...
        p5  p6];
    
else  %  the optimal column is from the second partiton of the Sylvester matrix poly g
    
    r_aboveG = mincol-n+t-2;
    r_belowG = m+n-2*t-mincol+2;
    
    p1 = zeros(r_aboveG,m+1);
    p2 = zeros(r_aboveG,n+1);
    
    p3 = zeros(n+1,m+1);
    
    G = RightDG_Build(n,m,t,mincol,theta);
    p5 = zeros(r_belowG,m+1);
    p6 = zeros(r_belowG,n+1);
    
    P = [p1     p2;...
        p3     G;...
        p5     p6];
    
end

end

function G = LeftDG_Build(m,n,t,mincol,theta)
% Build the matrix G, where G forms part of DP.
%
%
% Inputs
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% t : Degree of GCD 
%
% mincol :
%
% theta :
%
% Outputs.
%
% G : Matrix G


% Global Variables
global SETTINGS



switch SETTINGS.BOOL_LOG
    
    case 'y' % use log
        G = LeftDG_Build_log(m,n,t,mincol,theta);
        
    case 'n' % use nchoosek
        G = LeftDG_Build_nchoosek(m,n,t,mincol,theta);
        
    otherwise
        error('SETTINGS.BOOL_LOG must be either y or n')
end
end

function G = RightDG_Build(n,m,t,mincol,theta)
%
%
%
% Inputs.
%
% m : Degree of polynomial f
%
% n : Degree of polynomial g
%
% t : Degree of GCD
%
% mincol : index of column of S_{t}
%
% theta : Optimal value of theta


% Global Variables.
global SETTINGS

switch SETTINGS.BOOL_LOG
    case 'y' % use logs
        G = RightDG_Build_log(n,m,t,mincol,theta);
    case 'n' % use nchoosek
        warning('off','all');
        G = RightDG_Build_nchoosek(n,m,t,mincol,theta);
        warning('on','all');
    otherwise
        error('SETTINGS.BOOL_LOG must be either y or n')
end
end

function [G] = LeftDG_Build_nchoosek(m,n,t,num_mincolumn,theta)
%
%
% Inputs
%
% m : Degree of polynomial f
%
% n : Degree of polynomial g
%
% t : Degree of GCD
%
% mincol :
%
% theta :


% Global Variables
global SETTINGS

% Initalise vector G
G = zeros(1,m+1);

% Get column to be removed.
j = num_mincolumn-1;

for i = 0:1:m
    G(i+1) =  (theta.^(i)) ...
        .* nchoosek(i+j,j) ...
        .* nchoosek(m+n-t-(i+j),n-t-j);
    
end

% Create diagonal matrix G from values in vector G.
G = diag(G);


switch SETTINGS.BOOL_DENOM_SYL
    case 'y'
        % Included denominator in coefficient matrix.
        G = G./ nchoosek(m+n-t,n-t);
    case 'n'
        % Exclude
    otherwise
        error(' SETTINGS.BOOL_DENOM_SYL must be either y or n')
end

end

function [G] = LeftDG_Build_log(m,n,t,idxMinCol,theta)
% Get G \in\mathbb{R}^{(m+1)\times(m+1)}, where G is a diagonal matrix.
%
% Inputs
%
% m - Degree of polynomial f(x)
%
% n - Degree of polynomial g(x)
%
% t - Degree of GCD d(x)
%
% idxMinCol - Index of column removed from Sylvester Matrix
%
% theta - Optimal value of \theta

% Global Variables
global SETTINGS

% Initialise empty vector G.
G = zeros(1,m+1);

j = idxMinCol - 1;

for i = 0:1:m
    
    % Evaluate binomial coefficients in the numerator in logs
    Numerator_Eval_log = ...
        lnnchoosek(i+j,j) + ...
        lnnchoosek(m+n-t-(i+j),n-t-j);
    
    % Convert to standard form
    Numerator_Eval_exp = 10.^Numerator_Eval_log;
    
    %
    G(i+1) =  (theta.^(i)) .* Numerator_Eval_exp;
    
end

% Create diagonal matrix G from values in vector G.
G = diag(G);

switch SETTINGS.BOOL_DENOM_SYL
    case 'y'
        % Included denominator in coefficient matrix.
        
        % Evaluate binomial coefficient in denominator in logs
        Denom_Eval_log = lnnchoosek(m+n-t,n-t);
        
        % Conver to standard form
        Denom_Eval_exp = 10.^Denom_Eval_log;
        
        % Divide G by the denominator.
        G = G./Denom_Eval_exp;
    case 'n'
        % Exclude denominator from coefficient matrix
    otherwise
        error('SETTINGS.BOOL_DENOM_SYL must be either y or n')
        
end




end


function [G] = RightDG_Build_nchoosek(n,m,t,idxMinCol,theta)
%
%
% Inputs
%
% m - Degree of polynomial f(x)
%
% n - Degree of polynomial g(x)
%
% t - Degree of GCD
%
% idxMinCol -
%
% theta -



global SETTINGS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise empty vector G.
G = zeros(1,n+1);

j = idxMinCol-n+t-2;

for i=0:1:n
    
    G(i+1) = (theta.^(i)) .* nchoosek(i+j,j) .* nchoosek(m+n-t-(i+j),m-t-j) ;
    
end

% Create diagonal matrix G from vector G.
G = diag(G);

switch SETTINGS.BOOL_DENOM_SYL
    case 'y'
        % Denominator is included in coefficient matrix.
        G = G./ nchoosek(m+n-t,m-t);
    case 'n'
        % Denominator is excluded from the coefficient matrix
    otherwise
        error('SETTINGS.BOOL_DENOM_SYL must be either y or n')
        
end



end




function [G] = RightDG_Build_log(n,m,t,idxMinCol,theta)
%
%
% Inputs
%
% m - Degree of polynomial f
%
% n - Degree of polynomial g
%
% t - Degree of GCD
%
% mincol -
%
% theta -

% Global Variables
global SETTINGS

% Initialise the Right Matrix DG
G = zeros(1,n+1);


j = idxMinCol-n+t-2;


for i=0:1:n
    Numerator_Eval_log =...
        lnnchoosek(i+j,i) + ...
        lnnchoosek(m+n-t-(i+j),m-t-j);
    
    Numerator_Eval_exp = 10.^Numerator_Eval_log;
    
    G(i+1) = theta.^(i) .* Numerator_Eval_exp ;
    
end

% Create diagonal matrix G from vector G.
G = diag(G);

switch SETTINGS.BOOL_DENOM_SYL
    case 'y'
        % Include the denominator
        
        Denom_log = lnnchoosek(m+n-t,n);
        
        Denom_exp = 10.^Denom_log;
        
        G = G./ Denom_exp;
    case 'n'
        % Exclude the denominator
    otherwise
        error('SETTINGS.BOOL_DENOM_SYL must be either y or n')
end

end