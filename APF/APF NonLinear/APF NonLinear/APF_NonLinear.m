function [fx_lr, gx_lr, dx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = ...
    APF_NonLinear(fx,gx,ux,vx,i_alpha,i_theta,k)
% Refine Approximate Polynomial Factorisation by APF
%
%
% % Inputs
%
% fx : Coefficients of the polynomial f(x) in the Bernstein basis
%
% gx : Coefficients of the polynomial g(x) in the Bernstein basis
%
% ux : coefficients of the quotient polynomial u(x) in the Bernstein basis
%
% vx - Coefficients of the quotient polynomial v(x) in the Bernstein basis
%
% i_alpha - input \alpha
%
% i_theta - input \theta
%
% t : Calculated degree of GCD
%
% % Outputs
%
% fx_lr :
%
% gx_lr :
%
% dx_lr :
% 
% ux_lr :
%
% vx_lr :
%
% alpha_lr
%
% theta_lr



% Global Variables

global SETTINGS

% Initialise iteration index
ite = 1;


res_uw = zeros(1,1);
res_vw = zeros(1,1);
res_ux = zeros(1,1);
res_vx = zeros(1,1);
residual = zeros(1,1);

% Set initial values of alpha and theta
th(ite) = i_theta;
alpha(ite) = i_alpha;

% Get degree of polynomial f
m = GetDegree(fx);

% Get degree of polynomial g
n = GetDegree(gx);

% Initialise some useful vectors
vecm    = (0:1:m)';
vecn    = (0:1:n)';
veck    = (0:1:k)';

% Convert f and g to modified bernstein basis, excluding binomial
% coefficient
fw = GetWithThetas(fx,th);
gw = GetWithThetas(gx,th);

% Convert u and v to modified bernstein basis, excluding binomial
% coefficient
uw = GetWithThetas(ux,th);
vw = GetWithThetas(vx,th);

% Initialise S, such that sk = S * pt
S = (diag(th(ite).^vecm));

% Initialise T - Matrix such that tk = T * qt
T = (diag(th(ite).^vecn));

% Initialise zk - Structured perturbations of u and v
zk = zeros(m+n-2*k+2,1);

% Get partial derivatives of fw and gw with respect to theta
fw_wrt_theta = Differentiate_wrt_theta(fw,th(ite));
gw_wrt_theta = Differentiate_wrt_theta(gw,th(ite));



% Get partial derivatives of uw and vw with respect to theta.
Partial_uw_wrt_theta = Differentiate_wrt_theta(uw,th(ite));
Partial_vw_wrt_theta = Differentiate_wrt_theta(vw,th(ite));

% Get H^{-1} * C(u,v) * G
[HCG,H1C1G,H2C2G] = BuildHCG(uw,vw,k);




% Build HCG with respect to theta
[~,H1C1G_wrt_theta,H2C2G_wrt_theta] = ...
    BuildHCG(Partial_uw_wrt_theta,Partial_vw_wrt_theta,k);

%Build the RHS vector b = [fw,alpha.*gw]
bk = [fw ; alpha(ite).*gw];

dw = SolveAx_b(HCG,bk);
dx = GetWithoutThetas(dw,th(1));

% Get partial derivative of dw with respect to theta
dw_wrt_theta = Differentiate_wrt_theta(dw,th(ite));

% Get initial residual
res_vec = bk - ((HCG)*dw);

% Set some initial values
residual(ite)   = norm(res_vec);
zf               = zeros(m+1,1);
zg               = zeros(n+1,1);

% Values of LHS Perturbations
% Obtain structured perturbations sw of fw, and tw of gw
z_fw = GetWithThetas(zf,th);
z_gw = GetWithThetas(zg,th);

% Set initial values for the iterative process
z1_ux = zeros(length(uw),1);
z2_vx = zeros(length(vw),1);

% Construct the coefficient matrix in the equation that defines
% the constraint for the LSE problem.
HYk         = BuildHYQ_SNTLN(dx,m,n,th(ite));


% % Build the matrix C given by Hz Hp Hq Halpha Htheta1 Htheta2
H_z         = HYk;

H_p         = (-1)*S;

H_q         = -(alpha(ite))*T;

H_alpha     = -(gw + z_gw);

H_theta1    = -fw_wrt_theta+...
    (H1C1G_wrt_theta * dw)+...
    (H1C1G * dw_wrt_theta);

H_theta2    = (-(alpha(ite))*gw_wrt_theta)+...
    (H2C2G_wrt_theta * dw)+...
    (H2C2G * dw_wrt_theta);

C_temp      = ...
    [
    H_p,             zeros(m+1,n+1), zeros(m+1,1), H_theta1;...
    zeros(n+1,m+1),  H_q,            H_alpha,       H_theta2...
    ];

C       = [H_z , C_temp];

E       = eye(2*m+2*n-2*k+6);

fnew    = zeros(2*m+2*n-2*k+6,1);


ek = bk;

% Get the condition
condition(ite) = norm(res_vec)/norm(ek);

startpoint = [...
    zk;...
    zf;...
    zg;...
    alpha;...
    th];

yy = startpoint;

% Start the iterative procedure for the solution of the LSE problem.

while condition(ite) > (SETTINGS.MAX_ERROR_APF) && ite < SETTINGS.MAX_ITERATIONS_APF
    
    % Use the QR decomposition to solve the LSE problem.
    % min |y-p| subject to Cy=q
    y = LSE(E,fnew,C,res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % obtain the small changes.
    delta_zk        = y(1:m+n-2*k+2);
    delta_zf_k        = y(m+n-2*k+3:2*m+n-2*k+3);
    delta_zg_k        = y(2*m+n-2*k+4:2*m+2*n-2*k+4);
    delta_alpha     = y(end-1);
    delta_theta     = y(end);
    
    % Update variables zk, pk, qk, beta, theta
    zk          = zk + delta_zk;
    zf           = zf  + delta_zf_k;
    zg           = zg  + delta_zg_k;
    alpha(ite)  = alpha(ite-1) + delta_alpha;
    th(ite)     = th(ite-1) + delta_theta;
    
    % Update the iterative value of f and g
    fw = GetWithThetas(fx,th(ite));
    gw = GetWithThetas(gx,th(ite));
    
    % Update matrices S and T
    S = diag(th(ite).^vecm);
    T = diag(th(ite).^vecn);
    
    % Update the iterative value of GCD dw
    dw = GetWithThetas(dx,th(ite));
    
    % Update the iterative value of quotients uw and vw
    uw = GetWithThetas(ux,th(ite));
    vw = GetWithThetas(vx,th(ite));
    
    % Obtain partial derivatives of uw and vw with respect to theta
    uw_wrt_theta = Differentiate_wrt_theta(uw,th(ite));
    vw_wrt_theta = Differentiate_wrt_theta(vw,th(ite));
    
    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G
    [~,H1C1G,H2C2G] = BuildHCG(uw,vw,k);
    
    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G with respect to
    % theta
    [~,H1C1G_wrt_theta, H2C2G_wrt_theta] = BuildHCG(uw_wrt_theta,vw_wrt_theta,k);
    
    % Get perturbation vector, and seperate in to perturbations of f,
    % z1 and perturbations of g, z2
    z1_ux = zk(1:m-k+1);
    z2_vx = zk(m-k+2:m+n-2*k+2);
    
    % Obtain z1 and z2 in the modified bernstein basis.
    z_uw = GetWithThetas(z1_ux,th(ite));
    z_vw = GetWithThetas(z2_vx,th(ite));
    
    % obtain partial derivatives of z1 and z2 with respect to theta
    z2w_wrt_theta = Differentiate_wrt_theta(z_vw,th(ite));
    z1w_wrt_theta = Differentiate_wrt_theta(z_uw,th(ite));
    
    % Build Matrices H_{1}E_{1}(z1)G and H_{2}E_{2}(z2)G
    [~,H1E1G,H2E2G] = BuildHCG(z_uw,z_vw,k);
    
    % Calculate Partial derivatives of Matrices H_{1}E_{1}(z1)G and
    % H_{2}E_{2}(z2)G with respect to theta
    [~,H1E1G_wrt_theta,H2E2G_wrt_theta] = BuildHCG(z1w_wrt_theta,z2w_wrt_theta,k);
    
    % Obtain structured perturbations sw of fw, and tw of gw
    z_fw = GetWithThetas(zf,th(ite));
    z_gw = GetWithThetas(zg,th(ite));
    
    % Calculate partial derivatives of sw and tw with respect to theta
    s_wrt_theta = Differentiate_wrt_theta(z_fw,th(ite));
    t_wrt_theta = Differentiate_wrt_theta(z_gw,th(ite));
    
    % Calculate partial derivatives of fw and gw with respect to theta
    fw_wrt_theta = Differentiate_wrt_theta(fw,th(ite));
    gw_wrt_theta = Differentiate_wrt_theta(gx,th(ite));
    
    % Calculate partial derivative of dw with respect to theta
    dw_wrt_theta = Differentiate_wrt_theta(dw,th(ite));
    
    % Build Matrix C
    
    % Calculate H_z
    HYk = BuildHYQ_SNTLN(dx,m,n,th(ite));
    
    % Build H_z
    H_z = HYk;
    
    % Calculate H_p
    H_p = (-1)*S;
    
    % Calculate H_q
    H_q = -(alpha(ite))*T;
    
    % Calculate H_beta
    H_alpha = -(gw + z_gw);
    
    % Calculate H_theta1
    H_theta1 = -(fw_wrt_theta + s_wrt_theta)+...
        ((H1C1G_wrt_theta)*dw)+...
        ((H1E1G_wrt_theta)*dw)+...
        ((H1C1G)*dw_wrt_theta)+...
        ((H1E1G)*dw_wrt_theta);
    
    % Calculate H_theta2
    H_theta2 = -(alpha(ite))*(gw_wrt_theta + ...
        t_wrt_theta)+...
        ((H2C2G_wrt_theta)*dw)+...
        ((H2E2G_wrt_theta)*dw)+...
        ((H2C2G)*dw_wrt_theta)+...
        ((H2E2G)*dw_wrt_theta);
    
    C_temp = [H_p, zeros(m+1,n+1), zeros(m+1,1), H_theta1;
        zeros(n+1,m+1), H_q, H_alpha      , H_theta2];
    
    % Build Matrix C
    C = [H_z , C_temp];
    
    % Calculate Matrix H(C+E)G
    [HCEG,~,~] = BuildHCG(uw+z_uw,vw+z_vw,k);
    
    % Calculate the new residual
    res_vec = [fw + z_fw ;(alpha(ite)*(gw + z_gw))]-((HCEG)*dw);
    
    % Calculate the new right hand vector.
    ek = [fw + z_fw ;(alpha(ite)*(gw + z_gw))];
    
    
    % update vector of residual
    residual(ite) = norm(res_vec);
    
    % Update Condition scalar.
    condition(ite) = norm(res_vec)/norm(ek);
    
    % Update fnew
    fnew = -(yy - startpoint);
    
    % Edit 01/06/2015 16:30:00
    
    % Calculate the termination criterion in the modified
    % Bernstein basis..
    [res_uw(ite),res_vw(ite)] = Term_Criterion_APF(fw,gw,z_fw,z_gw,uw,vw,...
        dw,k,alpha(ite));
    
    % Repeat this calculation for the Bernstein basis. Transform the
    % variables from the modified Bernstein basis to the Bernstein
    % basis.
    
    fx_p = GetWithoutThetas(fw,th(ite));
    sx_p = GetWithoutThetas(z_fw,th(ite));
    gx_p = GetWithoutThetas(gw,th(ite));
    tx_p = GetWithoutThetas(z_gw,th(ite));
    ukx = ux + z1_ux;
    vkx = vx + z2_vx;
    dkx = dw./(th(ite).^veck);
    
    [res_ux(ite),res_vx(ite)] = Term_Criterion_APF(fx_p,gx_p,sx_p,...
        tx_p,ukx,vkx,dkx,k,1.0);
    
    
end





switch SETTINGS.PLOT_GRAPHS
    case 'y'
        % Plot the normalised residuals res_ux, res_vx, res_uw and res_vw.
        plotgraphs3(res_ux,res_vx,res_uw,res_vw);

        % Write out the number of iterations required and plot the values of
        % alpha, theta and the residual.
        plotgraphs4(alpha,th,residual);
    case 'n'
    otherwise
        error('err')
end

% Display number of iterations
LineBreakLarge();
fprintf('Iterations over approximate polynomial factorisation : %i \n', ite+1);
LineBreakLarge();

% Update values of quotients u and v,
ux_lr = ux + z1_ux;
vx_lr = vx + z2_vx;

% Update value of theta
theta_lr = th(ite);
alpha_lr = alpha(ite);

% Update value of common divisor dx
dx_lr = GetWithoutThetas(dw,th(ite));

% Edit 20/07/2015
fx_lr = GetWithoutThetas((fw + z_fw),th(ite));
gx_lr = GetWithoutThetas((gw + z_gw),th(ite));



switch SETTINGS.PLOT_GRAPHS
    case 'y'        
        figure_name = sprintf('%s : Residuals',mfilename);
        figure('name',figure_name)
        title('plotting residual')
        hold on
        plot(1:1:length(residual),(residual));
        
    case 'n'
    otherwise
        error('err')
end

end




