function [uw,vw ] = GetQuotients(Aki,cki,n,t,opt_col)
% obtain coefficients of quotient polynomials u and v

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs

% Aki :- LHS Sylvester matrix excluding optimal column

% cki :- RHS Column, removed from sylvester matrix S_{t}

% t :- Degree of the GCD.

% opt_col :- index of the removed column.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Outputs.

% uw and vw are in the modified bernstein form, u_{i}\theta^{i} and 
% v_{i}\theta^{i}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove the optimal column from the sylvester matrix to form Ak and ck,
% where Ak is the Sylvester Matrix with the column removed, and ck is the
% removed column.

% Display the condition number of the coefficient matrix DTQ (The Sylvester
% Matrix) used to calculate the solution vector x = [-v;u]

%display('Condition Number of DTQ used in finding quotient polynomials')
%disp(cond(Aki));
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





end