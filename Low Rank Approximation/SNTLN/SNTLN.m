function [ fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = ...
    SNTLN( fx, gx, i_alpha, i_th, k, idx_col)
% Given the polynomials f(x) and g(x) (normalised by geometric mean),
% obtain the low rank approximation of the Syvlester subresultant matrix 
% S(f,g). Return the perturbed polynomials f(x) and g(x)
%
% % Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x), in bernstein basis.
%
% gx : (Vector) Coefficients of polynomial g(x), in bernstein basis.
%
% i_alpha : (Float) Initial value of \alpha
%
% i_th :   (Float) Initial value of \theta
%
% t :   (Int) Degree of AGCD.
%
% idx_col : (Int) Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.
%
% % Outputs.
%
% fx_lr : (Vector) Coefficients of f(x) on output, in standard bernstein basis,
% including added structured perturbations.
%
% gx_lr : (Vector) Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.
%
% alpha_lr : (Float)
%
% theta_lr : (Float)
%
% ux_lr : (Vector) 
%
% vx_lr : (Vector)


global SETTINGS

% Set the initial iterations number
ite = 1;

% Set initial values of alpha and theta
theta(ite) = i_th;
alpha(ite) = i_alpha;

% Get degree of polynomials f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);
   
% Create the identity matrix I, the matrix M formed from I by removing the
% column equivalent to the optimal column for removal from the Sylvester 
% subresultant matrix, so that S_{t}(f,g)*M = A_{t}, where A_{t} is the
% Sylvester subresultant matrix with the column removed.
I = eye(m+n-2*k+2,m+n-2*k+2);
M = I;
M(:, idx_col) = [];

% Let e be the column removed from the identity matrix, such that
% S_{t}(f,g) * e gives the column c_{t}, where c_{t} is the optimal column
% removed from the Sylvester subresultant.
e = I(:, idx_col);

% Obtain polynomials in Modified Bernstein Basis, using initial values of
% alpha and theta.
fw = GetWithThetas(fx, theta(ite));
gw = GetWithThetas(gx, theta(ite));

% Form the Coefficient Matrix DTQ such that DTQ * x = [col]
DTQ = BuildDTQ(fw, alpha.*gw, k);

% Calculate the partial derivatives of fw and gw with respect to alpha
fw_wrt_alpha = zeros(m+1, 1);
alpha_gw_wrt_alpha = gw;

% Calculate the partial derivatives of fw and gw with respect to theta
fw_wrt_theta    = Differentiate_wrt_theta(fw, theta(ite));
gw_wrt_theta    = Differentiate_wrt_theta(gw, theta(ite));

% Calculate derivative of D_{k}T(f,g)Q_{k} with respect to alpha
DTQ_wrt_alpha = BuildDTQ(fw_wrt_alpha, alpha_gw_wrt_alpha, k);

% Calculate the derivative of D_{k}TQ_{k} with respect to theta
DTQ_wrt_theta = BuildDTQ(fw_wrt_theta, alpha(ite).*gw_wrt_theta, k);

% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.
zk = zeros(m+n+2, 1);
z_fx = zk(1 : m+1);
z_gx = zk(m+2 : end);

% Initilaise the derivative of D_{k}NQ_{k} wrt alpha.
DNQ_wrt_alpha   = zeros(m+n-k+1, m+n-(2*k)+2);

% Initilaise the derivative of D_{k}NQ_{k} wrt theta.
DNQ_wrt_theta   = zeros(m+n-k+1, m+n-(2*k)+2);

% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
hk_wrt_alpha     = DNQ_wrt_alpha*e;
hk_wrt_theta     = DNQ_wrt_theta*e;

%Calculate the matrix P.
DPQ = BuildDPG_SNTLN(m, n, k, alpha(ite), theta(ite), idx_col);

% Calculate the column of DTQ that is moved to the right hand side.
ck = (DTQ)*e;

% Calculate the remaining columns of the matrix.
At = (DTQ)*M;

% Calculate the derivatives wrt alpha and theta of the removed column.
ck_wrt_alpha     = DTQ_wrt_alpha *e;
ck_wrt_theta     = DTQ_wrt_theta *e;

% Calculate the initial estimate of x - the vector whcih contains the
% coefficients of the quotient polynomials u and v.
xk = SolveAx_b(At, ck);

% %
% Build Matrix Y

% Get vector x_{k}
x1 = xk(1 : idx_col-1) ;
x2 = xk(idx_col : end) ;
x = [x1; 0 ;x2] ;% Insert zero into vector

DYQ = BuildDYQ_SNTLN(x, m, n, k, alpha(ite), theta(ite));

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ck - (DTQ*M*xk);

% Set the intial value of E to the identity matrix
E = eye(2*m+2*n-2*k+5);

% Create the matrix D(T+N)Q
DTNQ = BuildDTQ(fw, alpha(ite).*gw, k);

% Create The matrix D(T+N)Q with respect to alpha
DTNQ_wrt_alpha = BuildDTQ(fw_wrt_alpha, alpha_gw_wrt_alpha, k);

% Create The matrix D(T+N)Q with respect to theta
DTNQ_wrt_theta = BuildDTQ(fw_wrt_theta, alpha(ite).*gw_wrt_theta, k);

% Create the matrix C for input into iteration

H_z     = DYQ - DPQ;

H_x     = DTNQ*M;

H_alpha  = DTNQ_wrt_alpha*M*xk - ...
    (ck_wrt_alpha + hk_wrt_alpha);

H_theta = DTNQ_wrt_theta*M*xk - ...
    (ck_wrt_theta + hk_wrt_theta);

C       = [H_z H_x H_alpha H_theta];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
        zk;...
        xk;...
        alpha(ite);...
        theta(ite)
    ];

yy              =   start_point;

% Set the initial value of vector p to be zero
f = zeros(2*m+2*n-2*k+5, 1);


% Set the termination criterion to a large value. It will be
% over written later.
condition(ite) = norm(res_vec)/norm(ck);

while condition(ite) >(SETTINGS.MAX_ERROR_SNTLN) &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
   
    
    % Use the QR decomposition to solve the LSE problem
    % min |y-p| subject to Cy=q
    y = LSE(E, f, C,res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % obtain the small changes
    delta_zk = y(1:m+n+2,1);
    delta_xk = y((m+n+3):(2*m+2*n-2*k+3),1);
    delta_alpha = y(2*m+2*n-2*k+4);
    delta_theta = y(2*m+2*n-2*k+5);
    
    % %
    % Update variables z_{k}, x_{k}, where z_{k} are perturbations in the
    % coefficients of f and g. x_{k} is the solution vector, containing
    % coefficients u and v.
    
    % Update the vector z
    zk = zk + delta_zk;
    
    % Update least squares solution
    xk = xk + delta_xk;
        
    % Update \alpha 
    alpha(ite) = alpha(ite-1) + delta_alpha;
    
    % Update \theta
    theta(ite) = theta(ite-1) + delta_theta;

    % Get f(\omega) from f(x) and g(\omega) from g(x)
    fw = GetWithThetas(fx,theta(ite));
    gw = GetWithThetas(gx,theta(ite));
    
    % Construct the subresultant matrix of DTQ.
    DTQ = BuildDTQ(fw,alpha(ite).*gw,k);
        
    % Get the partial derivative of f(\omega) with respect to \alpha
    fw_wrt_alpha    = zeros(m+1, 1);
    
    % Get the partial derivative of g(\omega) with respect to \alpha
    gw_wrt_alpha    = gw;
    
    % Get the partial derivative of f(\omega) with respect to \theta 
    % Get the partial derivative of g(\omega) with respect to \theta
    fw_wrt_theta    = Differentiate_wrt_theta(fw, theta(ite));
    gw_wrt_theta    = Differentiate_wrt_theta(gw, theta(ite));
    
    % Calculate the partial derivative of DTQ with respect to alpha.
    DTQ_wrt_alpha = BuildDTQ(fw_wrt_alpha, gw_wrt_alpha,k);
        
    % Calculate the partial derivative of DTQ with respect to theta.
    DTQ_wrt_theta = BuildDTQ(fw_wrt_theta,...
                alpha(ite).*gw_wrt_theta,k);
        
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ck = DTQ*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    ck_wrt_alpha     = DTQ_wrt_alpha * e;
    ck_wrt_theta     = DTQ_wrt_theta * e;
    
    % Create the vector of structured perturbations zf(x) and zg(x) applied
    % to F and G.
    z_fx      = zk(1 : m+1);
    z_gx      = zk(m+2 : end);
    
    % Get z_{f}(\omega) from z_{f}(x)
    z_fw = GetWithThetas(z_fx, theta(ite));
    
    % Get z_{g}(\omega) from z_{g}(x)
    z_gw = GetWithThetas(z_gx, theta(ite));
    
    % Calculate the derivatives of z_fw and z_gw with repect to alpha.
    zfw_wrt_alpha    = zeros(m+1, 1);
    zgw_wrt_alpha    = z_gw;
    
    % Calculate the derivatives of z_fw and z_gw with respect to theta.
    zfw_wrt_theta    = Differentiate_wrt_theta(z_fw, theta(ite));
    zgw_wrt_theta    = Differentiate_wrt_theta(z_gw, theta(ite));
    
    % Build the Coefficient Matrix DNQ, of structured perturbations, with
    % same structure as DTQ.
    DNQ = BuildDTQ(z_fw, alpha(ite).*z_gw, k);
    
    % Calculate the derivatives of DNQ with respect to alpha
    DNQ_wrt_alpha = BuildDTQ(zfw_wrt_alpha, zgw_wrt_alpha, k);
    
    
    % Calculate the derivatives of DNQ with respect to theta
    DNQ_wrt_theta = BuildDTQ(zfw_wrt_theta, alpha(ite).*zgw_wrt_theta, k);
    
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{t} the column of S_{t} moved to the RHS
    hk = DNQ * e;
    
    % Calculate the derivative of h_{t} with respect to alpha
    hk_alpha = DNQ_wrt_alpha * e;
    
    % Calculate the derivative of h_{t} with respect to theta
    hk_theta = DNQ_wrt_theta * e;
    
    % Build the matrix D_{t}(T+N)Q_{t}
    DTNQ = BuildDTQ(fw + z_fw, alpha(ite).*(gw + z_gw), k);
    
    % Calculate the paritial derivative of D_{t}(T+N)Q_{t} with respect to
    % alpha
    DTNQ_alpha = BuildDTQ(fw_wrt_alpha + zfw_wrt_alpha,...
                gw_wrt_alpha + zgw_wrt_alpha,k);
        
    % Calculate the paritial derivative of D_{k}(T+N)Q_{k} with respect to
    % theta
    DTNQ_theta = BuildDTQ(fw_wrt_theta + zfw_wrt_theta,...
                    alpha(ite).*(gw_wrt_theta + zgw_wrt_theta),k);
        
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    

    % %  
    % Build Matrix DY
    x1 = xk(1 : idx_col-1) ;
    x2 = xk(idx_col : end) ;
    x = [x1; 0 ;x2]; 
    
    % Build the matrix DY
    DYQ = BuildDYQ_SNTLN(x, m, n, k, alpha(ite), theta(ite));
    
    % Calculate the matrix DP where P is the matrix such that c = P[f;g]
    DPQ = BuildDPG_SNTLN(m, n, k, alpha(ite), theta(ite), idx_col);
    
    % Calculate the residual q and vector p.
    res_vec = (ck + hk) - (DTNQ * M * xk);
    
    % Create the matrix C. This is made up of four submatrices, HZ, Hx,
    % H_alpha and H_theta.
    
    Hz      = DYQ - DPQ;
    
    Hx      = DTNQ*M;
    
    H_alpha = DTNQ_alpha*M*xk - (ck_wrt_alpha + hk_alpha);
    
    H_theta = DTNQ_theta*M*xk - (ck_wrt_theta + hk_theta);
    
    C = [Hz,Hx,H_alpha,H_theta];  % the matrix C
       
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(res_vec)./ norm(ck + hk) ;
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    

    
end

% Plot graphs
Plot_STLN()
PlotThetas(theta)
PlotAlphas(alpha)

% Print the number of iterations required
LineBreakLarge()
fprintf([mfilename ' : ' sprintf('Iterations required for STLN : %i \n', ite)]);
LineBreakLarge()
SETTINGS.LOW_RANK_APPROX_REQ_ITE = ite;

% %
% Get polynomials for output

% Get the polynomial f(x) + \delta f(x)
fx_lr = fx + z_fx;

% Get the polynomial g(x) + \delta g(x)
gx_lr = gx + z_gx;

% % Get polynomials u(x) and v(x)

% Get u(x) and v(x) from x_ls
x = [xk(1:idx_col-1) ; -1 ; xk(idx_col:end)];

% Get the number of coefficients in v(x)
nCoefficients_vx = n-k+1;

% Get coefficients of v(\omega) from the vector x
vw_lr = x(1 : nCoefficients_vx);

% Get v(x) from v(\omega)
vx_lr = GetWithoutThetas(vw_lr, theta(ite));

% Get polynomial coefficients v(\omega) from the vector x
uw_lr = -1.*(x(nCoefficients_vx + 1:end));

% Get u(x) from u(\omega)
ux_lr = GetWithoutThetas(uw_lr,theta(ite));

% Get \alpha and \theta
alpha_lr = alpha(ite);
theta_lr = theta(ite);

fprintf('Change in theta \n')
display(theta(ite) - theta(1));
fprintf('Change in alpha \n')
display(alpha(ite) - alpha(1));

end



function PlotThetas(vTheta)
%
% % Inputs
%
% vTheta : vector of theta values for each iteration

figure_name = sprintf([mfilename ' : Theta variation over SNTLN']);
figure('name',figure_name)
hold on
plot(log10(vTheta),'-s','DisplayName','\theta')
xlabel('Iteration');
ylabel('log_{10}');
hold off

end

function PlotAlphas(vAlpha)


figure_name = sprintf([mfilename ' : Alpha variation over SNTLN']);
figure('name',figure_name)
hold on
plot(log10(vAlpha),'-s','DisplayName','\alpha')
xlabel('Iteration');
ylabel('log_{10}');
hold off

end








