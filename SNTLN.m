function [ fx_output,gx_output,alpha_output,theta_output,X_output] = ...
    SNTLN( fx_n,gx_n,initial_alpha,initial_theta,t,opt_col)
% Given the polynomials f(x) and g(x) (normalised by geometric mean),
% obtain the low rank approximation of the Syvlester matrix S(f,g) and the
% minimal perturbations such that S(f+\delta f, g+\delta g)x = 0.
%
%                             Inputs.
%
%
% fx_n :    Coefficients of polynomial f, in standard bernstein basis.
%
% gx_n :    Coefficients of polynomial g, in standard bernstein basis.
%
% initial_alpha :   Initial value of alpha
%
% initial_theta :   Initial value of theta
%
% t :   Degree of AGCD.
%
% opt_col : Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.
%
%
%
%                           Outputs.
%
%
% fx_output :- Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.
%
% gx_output :- Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.
%
% alpha_output :-
%
% theta_output :-
%
% X_output :-
%

%%
%                       Global Inputs

global MAX_ERROR_SNTLN
global MAX_ITERATIONS_SNTLN
global PLOT_GRAPHS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the initial iterations number
ite = 1;

% Set initial values of alpha and theta
theta(ite) = initial_theta;
alpha(ite) = initial_alpha;

% Get degree of polynomials f.
m = size(fx_n,1) - 1;

% Get degree of polynomial g.
n = size(gx_n,1) - 1;

% Generate some useful vectors
vecm = (0:1:m)';
vecn = (0:1:n)';

% Create the identity matrix I, the matrix M formed from I by removing the
% column equivalent to the optimal column for removal from the Sylvester 
% subresultant matrix, so that S_{t}(f,g)*M = A_{t}, where A_{t} is the
% Sylvester subresultant matrix with the column removed.
I = eye(m+n-2*t+2,m+n-2*t+2);

M = I;
M(:,opt_col) = [];

% Let e be the column removed from the identity matrix, such that
% S_{t}(f,g) * e gives the column c_{t}, where c_{t} is the optimal column
% removed from the Sylvester subresultant.
e = I(:,opt_col);

% Obtain polynomials in Modified Bernstein Basis, using initial values of
% alpha and theta.
fw = fx_n .* (theta(ite).^vecm);
gw = gx_n .* (theta(ite).^vecn);


% Form the Coefficient Matrix DTQ such that DTQ * x = [col]
DTQ = BuildDTQ(fw,gw,alpha,t);

% Calculate the partial derivatives of fw and gw with respect to alpha
Partial_fw_wrt_alpha            = zeros(m+1,1);
Partial_alpha_gw_wrt_alpha      = gw;

% Calculate the partial derivatives of fw and gw with respect to theta
Partial_fw_wrt_theta    = vecm.*fw./theta(ite);
Partial_gw_wrt_theta    = vecn.*gw./theta(ite);

% Calculate derivative of D_{k}T(f,g)Q_{k} with respect to alpha
Partial_DTQ_wrt_alpha = BuildDTQ(Partial_fw_wrt_alpha,Partial_alpha_gw_wrt_alpha,1,t);

% Calculate the derivative of D_{k}TQ_{k} with respect to theta
Partial_DTQ_wrt_theta = BuildDTQ(Partial_fw_wrt_theta,alpha(ite).*Partial_gw_wrt_theta,1,t);

% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.
zt = zeros(m+n+2,1);

% Initilaise the derivative of D_{k}NQ_{k} wrt alpha.
Partial_DNQ_wrt_alpha   = zeros(m+n-t+1,m+n-(2*t)+2);

% Initilaise the derivative of D_{k}NQ_{k} wrt theta.
Partial_DNQ_wrt_theta   = zeros(m+n-t+1,m+n-(2*t)+2);

% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
Partial_ht_wrt_alpha     = Partial_DNQ_wrt_alpha*e;
Partial_ht_wrt_theta     = Partial_DNQ_wrt_theta*e;

%Calculate the matrix P.
DP = BuildDP(m,n,theta(ite),opt_col,t);

% Calculate the column of DTQ that is moved to the right hand side.
ct = (DTQ)*e;

% Calculate the remaining columns of the matrix.
At = (DTQ)*M;

% Calculate the derivatives wrt alpha and theta of the removed column.
Partial_ct_wrt_alpha     = Partial_DTQ_wrt_alpha*e;
Partial_ct_wrt_theta     = Partial_DTQ_wrt_theta*e;

% Calculate the initial estimate of x - the vector whcih contains the
% coefficients of the quotient polynomials u and v.
x_ls = SolveAx_b(At,ct)

% Build Matrix Y
DY = BuildDY(m,n,t,opt_col,x_ls,alpha(ite),theta(ite));

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ct - (DTQ*M*x_ls);

% Set the initial value of vector p to be zero
p = zeros(2*m+2*n-2*t+5,1);

% Set the intial value of E to the identity matrix
E = eye(2*m+2*n-2*t+5);

% Create the matrix D(T+N)Q

DTNQ = BuildDTQ(fw,gw,alpha(ite),t);

% Create The matrix D(T+N)Q with respect to alpha
DTNQ_wrt_alpha = BuildDTQ(Partial_fw_wrt_alpha, Partial_alpha_gw_wrt_alpha,1,t);

% Create The matrix D(T+N)Q with respect to theta
DTNQ_wrt_theta = BuildDTQ(Partial_fw_wrt_theta, Partial_gw_wrt_theta,alpha(ite),t);

% Create the matrix C for input into iteration

H_z     = GetH_z(opt_col,n,t,DY,DP,alpha(ite));

H_x     = DTNQ*M;

H_alpha  = DTNQ_wrt_alpha*M*x_ls - ...
    (Partial_ct_wrt_alpha + Partial_ht_wrt_alpha);

H_theta = DTNQ_wrt_theta*M*x_ls - ...
    (Partial_ct_wrt_theta + Partial_ht_wrt_theta);

C       = [H_z H_x H_alpha H_theta];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
        zt;...
        x_ls;...
        alpha(ite);...
        theta(ite)
    ];

yy              =   start_point;

% Set the termination criterion to a large value. It will be
% over written later.
condition(ite) = norm(res_vec)/norm(ct);

while condition(ite) >(MAX_ERROR_SNTLN) &&  ite < MAX_ITERATIONS_SNTLN
   
    
    % Use the QR decomposition to solve the LSE problem
    % min |y-p| subject to Cy=q
    y = LSE(E,p,C,res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % obtain the small changes
    delta_zt        = y(1:m+n+2,1);
    delta_xt        = y((m+n+3):(2*m+2*n-2*t+3),1);
    delta_alpha     = y(2*m+2*n-2*t+4);
    delta_theta     = y(2*m+2*n-2*t+5);
    
    % Update variables z_{k}, x_{k}, where z_{k} are perturbations in the
    % coefficients of f and g. x_{k} is the solution vector, containing
    % coefficients u and v.
    zt   = zt   + delta_zt;
    x_ls = x_ls + delta_xt;
    
    % Update alpha and theta
    alpha(ite) = alpha(ite-1) + delta_alpha;
    theta(ite) = theta(ite-1) + delta_theta;
    
    % Obtain polynomials in modified bersntein basis a_{i}\theta^{i}
    fw = fx_n.*(theta(ite).^vecm);
    gw = gx_n.*(theta(ite).^vecn);
    
    % Construct the subresultant matrix of DTQ.
    DTQ = BuildDTQ(fw,gw,alpha(ite),t);
        
    % Calculate the partial derivatives of fw and gw with respect to alpha
    Partial_fw_wrt_alpha    = zeros(m+1,1);
    Partial_gw_wrt_alpha    = gw;
    
    % Calculate the partial derivatives of fw and gw with respect to theta
    Partial_fw_wrt_theta    = vecm.*fw./theta(ite);
    Partial_gw_wrt_theta    = vecn.*gw./theta(ite);
    
    % Calculate the Partial derivative of DTQ with respect to alpha.
    Partial_DTQ_wrt_alpha = BuildDTQ(Partial_fw_wrt_alpha, Partial_gw_wrt_alpha,1,t);
        
    % Calculate the partial derivative of DTQ with respect to theta.
    Partial_DTQ_wrt_theta = BuildDTQ(Partial_fw_wrt_theta,...
                Partial_gw_wrt_theta,alpha(ite),t);
        
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ct = DTQ*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    Partial_ct_wrt_alpha     = Partial_DTQ_wrt_alpha*e;
    Partial_ct_wrt_theta     = Partial_DTQ_wrt_theta*e;
    
    % Create the vector of structured perturbations zf(x) and zg(x) applied
    % to F and G.
    z_fx      = zt(1:m+1);
    z_gx      = zt(m+2:end);
    
    % Get zf(w) and zg(w)
    z_fw     = z_fx.*(theta(ite).^vecm);
    z_gw     = z_gx.*(theta(ite).^vecn);
    
    % Calculate the derivatives of z_fw and z_gw with repect to alpha.
    Partial_zfw_wrt_alpha    = zeros(m+1,1);
    Partial_zgw_wrt_alpha    = z_gw;
    
    % Calculate the derivatives of z_fw and z_gw with respect to theta.
    Partial_zfw_wrt_theta    = vecm.*(z_fw./theta(ite));
    Partial_zgw_wrt_theta    = vecn.*(z_gw./theta(ite));
    
    % Build the Coefficient Matrix DNQ, of structured perturbations, with
    % same structure as DTQ.
    DNQ = BuildDTQ(z_fw,z_gw,alpha(ite),t);
    
    % Calculate the derivatives of DNQ with respect to alpha
    Partial_DNQ_wrt_alpha = BuildDTQ(Partial_zfw_wrt_alpha,...
                Partial_zgw_wrt_alpha,1,t);
    
    
    % Calculate the derivatives of DNQ with respect to theta
    Partial_DNQ_wrt_theta = BuildDTQ(Partial_zfw_wrt_theta,...
                Partial_zgw_wrt_theta,alpha(ite),t);
    
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{t} the column of S_{t} moved to the RHS
    ht = DNQ*e;
    
    % Calculate the derivative of h_{t} with respect to alpha
    ht_alpha = Partial_DNQ_wrt_alpha*e;
    
    % Calculate the derivative of h_{t} with respect to theta
    ht_theta = Partial_DNQ_wrt_theta*e;
    
    % Build the matrix D_{t}(T+N)Q_{t}
    DTNQ = BuildDTQ(fw + z_fw, gw + z_gw ,alpha(ite),t);
    
    % Calculate the paritial derivative of D_{t}(T+N)Q_{t} with respect to
    % alpha
    DTNQ_alpha = BuildDTQ(Partial_fw_wrt_alpha + Partial_zfw_wrt_alpha,...
                Partial_gw_wrt_alpha + Partial_zgw_wrt_alpha,1,t);
        
    % Calculate the paritial derivative of D_{k}(T+N)Q_{k} with respect to
    % theta
    DTNQ_theta = BuildDTQ(Partial_fw_wrt_theta + Partial_zfw_wrt_theta,...
                (Partial_gw_wrt_theta + Partial_zgw_wrt_theta),alpha(ite),t);
        
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    DY = BuildDY(m,n,t,opt_col,x_ls,alpha(ite),theta(ite));
    
    % Calculate the matrix DP where P is the matrix such that c = P[f;g]
    DP = BuildDP(m,n,theta(ite),opt_col,t);
    
    % Calculate the residual q and vector p.
    res_vec = (ct+ht) - (DTNQ * M * x_ls);
    
    % Create the matrix C. This is made up of four submatrices, HZ, Hx,
    % H_alpha and H_theta.
    
    Hz      = GetH_z(opt_col,n,t,DY,DP,alpha(ite));
    
    Hx      = DTNQ*M;
    
    H_alpha = DTNQ_alpha*M*x_ls - (Partial_ct_wrt_alpha + ht_alpha);
    
    H_theta = DTNQ_theta*M*x_ls - (Partial_ct_wrt_theta + ht_theta);
    
    C = [Hz,Hx,H_alpha,H_theta];  % the matrix C
       
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(res_vec)./ norm(ct + ht) ;
    
    % Update fnew - used in LSE Problem.
    p = -(yy-start_point);
    

    
end

switch PLOT_GRAPHS
    case 'y'
        figure('name','SNTLN - Residuals')
        hold on
        title('residuals in SNTLN without constraints')
        xlabel('iterations')
        ylabel('residuals')
        plot(log10(condition),'-s')
        hold off
    case 'n'
    otherwise
        error('PLOT_GRAPHS must be either y or n');
end

if condition(ite) > condition(1)
    fprintf('SNTLN Failed to converge, default to input values\n')
    fx_output = fx_n;
    gx_output = gx_n;
    alpha_output = initial_alpha;
    theta_output = initial_theta;
    X_output = x_ls;
    return;
end

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.
fx_output = fx_n + zt(1:m+1);

gx_output = gx_n + zt(m+2:end);

X_output  = x_ls;

alpha_output = alpha(ite);

theta_output = theta(ite);



% Print the number of iterations
fprintf('--------------------------------------------------------------------------- \n')
fprintf('Iterations over Sylvester Matrix : %i \n', ite);
fprintf('--------------------------------------------------------------------------- \n')

end



function Hz = GetH_z(mincol,n,d,DY,DP,alpha)
if mincol<=(n-d+1)
    Hz=(DY-DP);
else
    Hz=DY-(alpha.*DP);
end
end













