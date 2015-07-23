function [PostAPF_fx, PostAPF_gx, PostAPF_dx, PostAPF_ux, PostAPF_vx, PostAPF_theta] = ...
    APF(fx,gx,ux,vx,initial_alpha,initial_theta,dx,t)
% Refine Approximate Polynomial Factorisation by newton raphson method iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                               Inputs

% fx - polynomial coefficients in bernstein basis

% gx - polynomial coefficients in bernstein basis

% ux - quotient polynomial coefficients

% vx - quotient polynomial coefficients

% alpha -

% theta -

% dx - initial approximation of GCD(fx,gx)

% t - Calculated degree of GCD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Global Variables

% bool_plotgraphs (bool)
% 0 - Dont plot graphs
% 1 - Plot graphs of LSE problem, alpha and theta
global bool_plotgraphs

% set limit of error in LSE problem
global max_error


% set maximum number of iterations for LSE problem
global max_iterations

% BuildMethod (bool)
%   0 - Naive build method D,T,Q calculated independently
%   1 - New build method, build elementwise
global Bool_APFBuildMethod

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise iteration index
ite = 1;


res_uw = zeros(1,1);
res_vw = zeros(1,1);
res_ux = zeros(1,1);
res_vx = zeros(1,1);
residual = zeros(1,1);

% Set initial values of alpha and theta
theta(ite) = initial_theta;
alpha(ite) = initial_alpha;

% Get degree of polynomial f
m = length(fx) - 1;

% Get degree of polynomial g
n = length(gx) - 1;

% Initialise some useful vectors
vecm    = (0:1:m)';
vecn    = (0:1:n)';
vecmk   = (0:1:m-t)';
vecnk   = (0:1:n-t)';
veck    = (0:1:t)';


% Convert f and g to modified bernstein basis, excluding binomial
% coefficient
fw = fx.*(theta(ite).^vecm);
gw = gx.*(theta(ite).^vecn);

% Convert u and v to modified bernstein basis, excluding binomial
% coefficient
uw = ux.*(theta(ite).^vecmk);
vw = vx.*(theta(ite).^vecnk);

% Convert d to modified bernstein basis, excluding binomial coefficient
dw = dx.*(theta(ite).^veck);

% Initialise S, such that sk = S * pt
S = (diag(theta(ite).^vecm));

% Initialise T - Matrix such that tk = T * qt
T = (diag(theta(ite).^vecn));

% Initialise zk - Structured perturbations of u and v
zk = zeros(m+n-2*t+2,1);

% Get partial derivatives of fw and gw with respect to theta
Partial_fw_wrt_theta = vecm .*fw ./ theta(ite);
Partial_gw_wrt_theta = vecn .*gw ./ theta(ite);

% Get partial derivative of dw with respect to theta
Partial_dw_wrt_theta     = veck.*dw./theta(ite);

% Get partial derivatives of uw and vw with respect to theta.
Partial_uw_wrt_theta = vecmk .*uw./theta(ite);
Partial_vw_wrt_theta = vecnk .*vw./theta(ite);

% Build H_{t}C(f,g)G_{t}
switch Bool_APFBuildMethod
    case 0
        H1 = BuildH1(m);
        H2 = BuildH1(n);
        C1 = BuildC1Partition(uw,t);
        C2 = BuildC1Partition(vw,t);
        G = BuildG(t);
        
        HCG = blkdiag(H1,H2)*[C1;C2]*G;
        H1C1G = H1*C1*G;
        H2C2G = H2*C2*G;
    case 1
        H1C1G = BuildHCGPart(uw,t);
        H2C2G = BuildHCGPart(vw,t);
        HCG = [H1C1G ; H2C2G ];
end


% Build HCG with respect to theta
switch Bool_APFBuildMethod
    case 0
        C1_wrt_theta = BuildC1Partition(Partial_uw_wrt_theta,t);
        C2_wrt_theta = BuildC1Partition(Partial_vw_wrt_theta,t);
        
        H1C1G_wrt_theta = H1*C1_wrt_theta*G;
        H2C2G_wrt_theta = H2*C2_wrt_theta*G;
        
    case 1
        H1C1G_wrt_theta = BuildHCGPart(Partial_uw_wrt_theta,t);
        H2C2G_wrt_theta = BuildHCGPart(Partial_vw_wrt_theta,t);
        
end

%Build the RHS vector b = [fw,alpha.*gw]
bk = [fw ; alpha(ite).*gw];

% Get initial residual
rk = bk-((HCG)*dw);

% Set some initial values
gnew            = rk;
residual(ite)   = norm(gnew);
p               = zeros(m+1,1);
q               = zeros(n+1,1);

% Values of LHS Perturbations
% Obtain structured perturbations sw of fw, and tw of gw
sw   = p.*(theta(ite).^vecm);
tw   = (q.*(theta(ite).^vecn));


% Set initial values for the iterative process
z1x = zeros(length(uw),1);
z2x = zeros(length(vw),1);



% Construct the coefficient matrix in the equation that defines
% the constraint for the LSE problem.
HYk         = BuildHYQ(dx,m,n,alpha(ite),theta(ite));

% Update the iterative value of f and g


H_z         = HYk;
H_p         = (-1)*S;
H_q         = -(alpha(ite))*T;
H_alpha     = -(gw + tw);

H_theta1    = -Partial_fw_wrt_theta+...
    (H1C1G_wrt_theta * dw)+...
    (H1C1G * Partial_dw_wrt_theta);

H_theta2    = (-(alpha(ite))*Partial_gw_wrt_theta)+...
    (H2C2G_wrt_theta * dw)+...
    (H2C2G * Partial_dw_wrt_theta);

C_temp      = [H_p,             zeros(m+1,n+1), zeros(m+1,1), H_theta1;...
               zeros(n+1,m+1),  H_q,            H_alpha,       H_theta2];

C       = [H_z , C_temp];

E       = eye(2*m+2*n-2*t+6);

fnew    = zeros(2*m+2*n-2*t+6,1);


ek = bk;

condition = norm(rk)/norm(ek);

startpoint = [zk;p;q;alpha;theta];

yy = startpoint;

% Start the iterative procedure for the solution of the LSE problem.

while condition > (max_error) && ite < max_iterations 
    
    % Use the QR decomposition to solve the LSE problem.
    % min |y-p| subject to Cy=q
    y = LSE(E,fnew,C,gnew);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % obtain the small changes.
    delta_zk        = y(1:m+n-2*t+2);
    delta_pk        = y(m+n-2*t+3:2*m+n-2*t+3);
    delta_qk        = y(2*m+n-2*t+4:2*m+2*n-2*t+4);
    delta_alpha     = y(end-1);
    delta_theta     = y(end);

    
    % Update variables zk, pk, qk, beta, theta
    zk          = zk + delta_zk;
    p           = p  + delta_pk;
    q           = q  + delta_qk;
    alpha(ite)  = alpha(ite-1) + delta_alpha;
    theta(ite)  = theta(ite-1) + delta_theta;
    
    % Update the iterative value of f and g
    fw = fx.*(theta(ite).^vecm);
    gw = gx.*(theta(ite).^vecn);
    
    % Update matrices S and T
    S = diag(theta(ite).^vecm);
    T = diag(theta(ite).^vecn);
    
    % Update vectors of thetas
    ok1 = theta(ite).^vecmk;
    ok2 = theta(ite).^vecnk;
    ok3 = theta(ite).^veck;
    
    % Update the iterative value of GCD dw
    dw = dx.*ok3;
    
    % Update the iterative value of quotients uw and vw
    uw = ux.*ok1;
    vw = vx.*ok2;
    
    % Obtain partial derivatives of uw and vw with respect to theta
    uw_wrt_theta = vecmk .* uw ./ theta(ite);
    vw_wrt_theta = vecnk .* vw ./ theta(ite);
    
    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G
    switch Bool_APFBuildMethod
        case 0
            C1 = BuildC1Partition(uw,t);
            C2 = BuildC1Partition(vw,t);
            H1C1G = H1*C1*G;
            H2C2G = H2*C2*G;
        case 1
            H1C1G = BuildHCGPart(uw,t);
            H2C2G = BuildHCGPart(vw,t);
    end
    
    
    
    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G with respect to
    % theta
    switch Bool_APFBuildMethod
        case 0
            C1_wrt_theta = BuildC1Partition(uw_wrt_theta,t);
            C2_wrt_theta = BuildC1Partition(vw_wrt_theta,t);
            
            H1C1G_wrt_theta = H1*C1_wrt_theta*G;
            H2C2G_wrt_theta = H2*C2_wrt_theta*G;
        case 1
            
            H1C1G_wrt_theta = BuildHCGPart(uw_wrt_theta,t);
            H2C2G_wrt_theta = BuildHCGPart(vw_wrt_theta,t);
            
    end
    
    % Get perturbation vector, and seperate in to perturbations of f,
    % z1 and perturbations of g, z2
    z1x = zk(1:m-t+1);
    z2x = zk(m-t+2:m+n-2*t+2);
    
    % Obtain z1 and z2 in the modified bernstein basis.
    z1w = z1x .* ok1;
    z2w = z2x .* ok2;
    
    % obtain partial derivatives of z1 and z2 with respect to theta
    z2w_wrt_theta = vecnk .* z2w ./ theta(ite);
    z1w_wrt_theta = vecmk .* z1w ./ theta(ite);
    
    % Build Matrices H_{1}E_{1}(z1)G and H_{2}E_{2}(z2)G
    switch Bool_APFBuildMethod
        case 0
            E1 = BuildC1Partition(z1w,t);
            E2 = BuildC1Partition(z2w,t);
            H1E1G = H1*E1*G;
            H2E2G = H2*E2*G;
            
        case 1
            
            H1E1G = BuildHCGPart(z1w,t);
            H2E2G = BuildHCGPart(z2w,t);
            
    end
    
    % Calculate Partial derivatives of Matrices H_{1}E_{1}(z1)G and
    % H_{2}E_{2}(z2)G with respect to theta
    
    switch Bool_APFBuildMethod
        case 0
            E1_wrt_theta = BuildC1Partition(z1w_wrt_theta,t);
            E2_wrt_theta = BuildC1Partition(z2w_wrt_theta,t);
            
            H1E1G_wrt_theta = H1 * E1_wrt_theta * G;
            H2E2G_wrt_theta = H2 * E2_wrt_theta * G;
            
        case 1
            H1E1G_wrt_theta = BuildHCGPart(z1w_wrt_theta,t);
            H2E2G_wrt_theta = BuildHCGPart(z2w_wrt_theta,t);
    end
    
    
    % Obtain structured perturbations sw of fw, and tw of gw
    sw   = p.*(theta(ite).^vecm);
    tw   = q.*(theta(ite).^vecn);
    
    % Calculate partial derivatives of sw and tw with respect to theta
    Partial_s_wrt_theta = vecm.*(p.*(theta(ite).^(vecm-1)));
    Partial_t_wrt_theta = vecn.*(q.*(theta(ite).^(vecn-1)));
    
    % Calculate partial derivatives of fw and gw with respect to theta
    Partial_fw_wrt_theta = vecm.*(fx.*(theta(ite).^(vecm-1)));
    Partial_gw_wrt_theta = vecn.*(gx.*(theta(ite).^(vecn-1)));
    
    % Calculate partial derivative of dw with respect to theta
    Partial_dw_wrt_theta = veck.*(dx.*(theta(ite).^(veck-1)));
    
    % Calculate H_z
    HYk = BuildHYQ(dx,m,n,alpha(ite),theta(ite));
    
    H_z = HYk;
    
    % Calculate H_p
    H_p = (-1)*S;
    
    % Calculate H_q
    H_q = -(alpha(ite))*T;
    
    % Calculate H_beta
    H_alpha = -(gw + tw);
    
    % Calculate H_theta1
    H_theta1 = -(Partial_fw_wrt_theta + Partial_s_wrt_theta)+...
        ((H1C1G_wrt_theta)*dw)+...
        ((H1E1G_wrt_theta)*dw)+...
        ((H1C1G)*Partial_dw_wrt_theta)+...
        ((H1E1G)*Partial_dw_wrt_theta);
    % Calculate H_theta2
    H_theta2 = -(alpha(ite))*(Partial_gw_wrt_theta + ...
        Partial_t_wrt_theta)+...
        ((H2C2G_wrt_theta)*dw)+...
        ((H2E2G_wrt_theta)*dw)+...
        ((H2C2G)*Partial_dw_wrt_theta)+...
        ((H2E2G)*Partial_dw_wrt_theta);
    
    C_temp = [H_p, zeros(m+1,n+1), zeros(m+1,1), H_theta1;
        zeros(n+1,m+1), H_q, H_alpha      , H_theta2];
    
    % Build Matrix C
    C = [H_z , C_temp];
    
    % Calculate Matrix H(C+E)G
    switch Bool_APFBuildMethod
        case 0
            H1C1E1G = H1 * (C1+E1) * G;
            H2C2E2G = H2 * (C2+E2) * G;
            HCEG = [H1C1E1G; H2C2E2G];
        case 1
            H1C1E1G = BuildHCGPart(uw+z1w,t);
            H2C2E2G = BuildHCGPart(vw+z2w,t);
            HCEG    = [H1C1E1G; H2C2E2G];
    end
    % Calculate the new residual
    rk = [fw + sw ;(alpha(ite)*(gw + tw))]-((HCEG)*dw);
    
    % Calculate the new right hand vector.
    ek = [fw + sw ;(alpha(ite)*(gw + tw))];
    
    
    % Update gnew - used in LSE Problem
    gnew = rk;
    
    % update vector of residual
    residual(ite) = norm(gnew);
    
    % Update Condition scalar.
    condition = norm(rk)/norm(ek);
    
    % Update fnew
    fnew = -(yy - startpoint);
    
    % Edit 01/06/2015 16:30:00
    
    % Calculate the termination criterion in the modified
    % Bernstein basis..
    [res_uw(ite),res_vw(ite)] = Term_Criterion_APF(fw,gw,sw,tw,uw,vw,...
        dw,t,alpha(ite));
    
    % Repeat this calculation for the Bernstein basis. Transform the
    % variables from the modified Bernstein basis to the Bernstein
    % basis.
    fx_p = fw./(theta(ite).^vecm); 
    sx_p = sw./(theta(ite).^vecm); 
    gx_p = gw./(theta(ite).^vecn);
    tx_p = tw./(theta(ite).^vecn);
    ukx = ux + z1x;   
    vkx = vx + z2x;
    dkx = dw./(theta(ite).^veck);
    
    [res_ux(ite),res_vx(ite)] = Term_Criterion_APF(fx_p,gx_p,sx_p,...
        tx_p,ukx,vkx,dkx,t,1.0);
    % Edit End
    
end


if ite == max_iterations 
    fprintf('APF failed to converge after %i iterations\n',ite)
    PostAPF_dx = dx; 
    PostAPF_ux = ux;
    PostAPF_vx = vx;
    PostAPF_theta = initial_theta;
    
    % Edit 20/07/2015
    PostAPF_fx = fx;
    PostAPF_gx = gx;
    return 
    
end

try
% EDIT 01/06/2015 16:34:00
switch bool_plotgraphs
    case 1
% Plot the normalised residuals res_ux, res_vx, res_uw and res_vw.
plotgraphs3(res_ux,res_vx,res_uw,res_vw);

% Write out the number of iterations required and plot the values of
% alpha, theta and the residual.
plotgraphs4(alpha,theta,residual);
end
catch
    fprintf('Can not perform plotgraphs3 and plotgraphs4, not enough iterations performed\n')
end
% END EDIT 01/06/2015 16:34:00


% Display number of iterations
fprintf('--------------------------------------------------------------------------- \n')
fprintf('Iterations over approximate polynomial factorisation : %i \n', ite+1);
fprintf('--------------------------------------------------------------------------- \n')

% Update values of quotients u and v,
PostAPF_ux = ux + z1x;
PostAPF_vx = vx + z2x;

% Update value of theta
PostAPF_theta = theta(ite);

% Update value of common divisor dx
PostAPF_dx = dw./(theta(ite).^veck);

% Edit 20/07/2015
PostAPF_fx = (fw + sw) ./ (theta(ite).^vecm);
PostAPF_gx = (gw + tw) ./ (theta(ite).^vecn);



switch bool_plotgraphs
    case 1
        
        plot_residuals = 1;
        plot_thetas = 1;
        plot_betas = 1;
        
        switch plot_residuals
            case 1
                figure()
                title('plotting residual')
                hold on
                plot(1:1:length(residual),(residual));
        end
        switch plot_thetas
            case 1
                figure()
                title('plotting theta')
                hold on
                plot(1:1:length(theta),(theta));
        end
        switch plot_betas
            case 1
                figure()
                title('plotting beta');
                hold on
                plot(1:1:length(alpha),(alpha));
        end
end

end




