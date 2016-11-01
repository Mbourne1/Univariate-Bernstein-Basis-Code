function [ fx_output,gx_output,alpha_output,theta_output,X_output] = ...
    SNTLN( fx_n,gx_n,i_alpha,i_th,t,opt_col)
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
% i_alpha :   Initial value of alpha
%
% i_th :   Initial value of theta
%
% t :   Degree of AGCD.
%
% opt_col : Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.
%
%
%
% % Outputs.
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

% %
% Global Inputs

global SETTINGS



% Set the initial iterations number
ite = 1;

% Set initial values of alpha and theta
th(ite) = i_th;
alpha(ite) = i_alpha;

% Get degree of polynomials f.
m = GetDegree(fx_n);

% Get degree of polynomial g.
n = GetDegree(gx_n);

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
fw = GetWithThetas(fx_n,th(ite));
gw = GetWithThetas(gx_n,th(ite));

% Form the Coefficient Matrix DTQ such that DTQ * x = [col]
DTQ = BuildDTQ(fw,alpha.*gw,t);

% Calculate the partial derivatives of fw and gw with respect to alpha
fw_wrt_alpha            = zeros(m+1,1);
alpha_gw_wrt_alpha      = gw;

% Calculate the partial derivatives of fw and gw with respect to theta
fw_wrt_theta    = Differentiate_wrt_theta(fw,th(ite));
gw_wrt_theta    = Differentiate_wrt_theta(gw,th(ite));

% Calculate derivative of D_{k}T(f,g)Q_{k} with respect to alpha
DTQ_wrt_alpha = BuildDTQ(fw_wrt_alpha,alpha_gw_wrt_alpha,t);

% Calculate the derivative of D_{k}TQ_{k} with respect to theta
DTQ_wrt_theta = BuildDTQ(fw_wrt_theta,alpha(ite).*gw_wrt_theta,t);

% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.
zt = zeros(m+n+2,1);

% Initilaise the derivative of D_{k}NQ_{k} wrt alpha.
DNQ_wrt_alpha   = zeros(m+n-t+1,m+n-(2*t)+2);

% Initilaise the derivative of D_{k}NQ_{k} wrt theta.
DNQ_wrt_theta   = zeros(m+n-t+1,m+n-(2*t)+2);

% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
ht_wrt_alpha     = DNQ_wrt_alpha*e;
ht_wrt_theta     = DNQ_wrt_theta*e;

%Calculate the matrix P.
DP = BuildDP(m,n,th(ite),opt_col,t);

% Calculate the column of DTQ that is moved to the right hand side.
ct = (DTQ)*e;

% Calculate the remaining columns of the matrix.
At = (DTQ)*M;

% Calculate the derivatives wrt alpha and theta of the removed column.
ct_wrt_alpha     = DTQ_wrt_alpha*e;
ct_wrt_theta     = DTQ_wrt_theta*e;

% Calculate the initial estimate of x - the vector whcih contains the
% coefficients of the quotient polynomials u and v.
x_ls = SolveAx_b(At,ct);

% Build Matrix Y
xa = x_ls(1:opt_col-1) ;
xb = x_ls(opt_col:end) ;
x = [xa; 0 ;xb] ;% Insert zero into vector

% Get number of cols in left partition
nCols_leftPart = (n-t+1);

% Get number of cols in right partition
nCols_rightPart = (m-t+1);

% Get the coefficients of xv and xu
xv = x(1:nCols_leftPart);
xu = x(nCols_leftPart+1:nCols_leftPart+nCols_rightPart);


DY = BuildDY(xu,xv,t,alpha(ite),th(ite));

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ct - (DTQ*M*x_ls);


% Set the intial value of E to the identity matrix
E = eye(2*m+2*n-2*t+5);

% Create the matrix D(T+N)Q
DTNQ = BuildDTQ(fw,alpha(ite).*gw,t);

% Create The matrix D(T+N)Q with respect to alpha
DTNQ_wrt_alpha = BuildDTQ(fw_wrt_alpha, alpha_gw_wrt_alpha,t);

% Create The matrix D(T+N)Q with respect to theta
DTNQ_wrt_theta = BuildDTQ(fw_wrt_theta, alpha(ite).*gw_wrt_theta,t);

% Create the matrix C for input into iteration

H_z     = GetH_z(opt_col,n,t,DY,DP,alpha(ite));

H_x     = DTNQ*M;

H_alpha  = DTNQ_wrt_alpha*M*x_ls - ...
    (ct_wrt_alpha + ht_wrt_alpha);

H_theta = DTNQ_wrt_theta*M*x_ls - ...
    (ct_wrt_theta + ht_wrt_theta);

C       = [H_z H_x H_alpha H_theta];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
        zt;...
        x_ls;...
        alpha(ite);...
        th(ite)
    ];

yy              =   start_point;

% Set the initial value of vector p to be zero
f = zeros(2*m+2*n-2*t+5,1);


% Set the termination criterion to a large value. It will be
% over written later.
condition(ite) = norm(res_vec)/norm(ct);

while condition(ite) >(SETTINGS.MAX_ERROR_SNTLN) &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
   
    
    % Use the QR decomposition to solve the LSE problem
    % min |y-p| subject to Cy=q
    y = LSE(E,f,C,res_vec);
    
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
    th(ite) = th(ite-1) + delta_theta;
    
    % Obtain polynomials in modified bersntein basis a_{i}\theta^{i}
    fw = GetWithThetas(fx_n,th(ite));
    gw = GetWithThetas(gx_n,th(ite));
    
    % Construct the subresultant matrix of DTQ.
    DTQ = BuildDTQ(fw,alpha(ite).*gw,t);
        
    % Calculate the partial derivatives of fw and gw with respect to alpha
    fw_wrt_alpha    = zeros(m+1,1);
    gw_wrt_alpha    = gw;
    
    % Calculate the partial derivatives of fw and gw with respect to theta
    fw_wrt_theta    = Differentiate_wrt_theta(fw,th(ite));
    gw_wrt_theta    = Differentiate_wrt_theta(gw,th(ite));
    
    % Calculate the Partial derivative of DTQ with respect to alpha.
    DTQ_wrt_alpha = BuildDTQ(fw_wrt_alpha, gw_wrt_alpha,t);
        
    % Calculate the partial derivative of DTQ with respect to theta.
    DTQ_wrt_theta = BuildDTQ(fw_wrt_theta,...
                alpha(ite).*gw_wrt_theta,t);
        
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ct = DTQ*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    ct_wrt_alpha     = DTQ_wrt_alpha*e;
    ct_wrt_theta     = DTQ_wrt_theta*e;
    
    % Create the vector of structured perturbations zf(x) and zg(x) applied
    % to F and G.
    z_fx      = zt(1:m+1);
    z_gx      = zt(m+2:end);
    
    % Get zf(w) and zg(w)
    z_fw = GetWithThetas(z_fx,th(ite));
    z_gw = GetWithThetas(z_gx,th(ite));
    
    % Calculate the derivatives of z_fw and z_gw with repect to alpha.
    zfw_wrt_alpha    = zeros(m+1,1);
    zgw_wrt_alpha    = z_gw;
    
    % Calculate the derivatives of z_fw and z_gw with respect to theta.
    zfw_wrt_theta    = Differentiate_wrt_theta(z_fw,th(ite));
    zgw_wrt_theta    = Differentiate_wrt_theta(z_gw,th(ite));
    
    % Build the Coefficient Matrix DNQ, of structured perturbations, with
    % same structure as DTQ.
    DNQ = BuildDTQ(z_fw,alpha(ite).*z_gw,t);
    
    % Calculate the derivatives of DNQ with respect to alpha
    DNQ_wrt_alpha = BuildDTQ(zfw_wrt_alpha,...
                zgw_wrt_alpha,t);
    
    
    % Calculate the derivatives of DNQ with respect to theta
    DNQ_wrt_theta = BuildDTQ(zfw_wrt_theta,...
                alpha(ite).*zgw_wrt_theta,t);
    
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{t} the column of S_{t} moved to the RHS
    ht = DNQ*e;
    
    % Calculate the derivative of h_{t} with respect to alpha
    ht_alpha = DNQ_wrt_alpha*e;
    
    % Calculate the derivative of h_{t} with respect to theta
    ht_theta = DNQ_wrt_theta*e;
    
    % Build the matrix D_{t}(T+N)Q_{t}
    DTNQ = BuildDTQ(fw + z_fw, alpha(ite).*(gw + z_gw),t);
    
    % Calculate the paritial derivative of D_{t}(T+N)Q_{t} with respect to
    % alpha
    DTNQ_alpha = BuildDTQ(fw_wrt_alpha + zfw_wrt_alpha,...
                gw_wrt_alpha + zgw_wrt_alpha,t);
        
    % Calculate the paritial derivative of D_{k}(T+N)Q_{k} with respect to
    % theta
    DTNQ_theta = BuildDTQ(fw_wrt_theta + zfw_wrt_theta,...
                    alpha(ite).*(gw_wrt_theta + zgw_wrt_theta),t);
        
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    %x_ls = SolveAx_b(DTNQ * M ,ct+ht);

    % Build Matrix Y
    xa = x_ls(1:opt_col-1) ;
    xb = x_ls(opt_col:end) ;
    x = [xa; 0 ;xb] ;% Insert zero into vector
    
    % Get number of cols in left partition
    nCols_leftPart = (n-t+1);
    
    % Get number of cols in right partition
    nCols_rightPart = (m-t+1);
    
    % Get the coefficients of xv and xu
    xv = x(1:nCols_leftPart);
    xu = x(nCols_leftPart+1:nCols_leftPart+nCols_rightPart);
    
    % Build the matrix DY
    DY = BuildDY(xu,xv,t,alpha(ite),th(ite));
    
    % Calculate the matrix DP where P is the matrix such that c = P[f;g]
    DP = BuildDP(m,n,th(ite),opt_col,t);
    
    % Calculate the residual q and vector p.
    res_vec = (ct+ht) - (DTNQ * M * x_ls);
    
    % Create the matrix C. This is made up of four submatrices, HZ, Hx,
    % H_alpha and H_theta.
    
    Hz      = GetH_z(opt_col,n,t,DY,DP,alpha(ite));
    
    Hx      = DTNQ*M;
    
    H_alpha = DTNQ_alpha*M*x_ls - (ct_wrt_alpha + ht_alpha);
    
    H_theta = DTNQ_theta*M*x_ls - (ct_wrt_theta + ht_theta);
    
    C = [Hz,Hx,H_alpha,H_theta];  % the matrix C
       
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(res_vec)./ norm(ct + ht) ;
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    

    
end


switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf('%s : Residuals',mfilename);
        figure('name',figure_name)
        hold on
        title('Residuals in SNTLN without constraints')
        xlabel('Iterations')
        ylabel('log_{10} Residuals')
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
    alpha_output = i_alpha;
    theta_output = i_th;
    X_output = x_ls;
    return;
end

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.
fx_output = fx_n + z_fx;

gx_output = gx_n + z_gx;

X_output  = x_ls;

alpha_output = alpha(ite);

theta_output = th(ite);



% Print the number of iterations

fprintf([mfilename ' : ' sprintf('Iterations required for STLN : %i \n', ite)]);

end



function Hz = GetH_z(mincol,n,d,DY,DP,alpha)
if mincol<=(n-d+1)
    Hz=(DY-DP);
else
    Hz=DY-(alpha.*DP);
end
end













