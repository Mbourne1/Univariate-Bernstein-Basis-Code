function [PostAPF_fx, PostAPF_gx, PostAPF_dx, PostAPF_ux, PostAPF_vx, PostAPF_theta] = ...
    APF(fx,gx,ux,vx,initial_alpha,initial_theta,dx,t)

% Refine Approximate Polynomial Factorisation by

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               Inputs
%
% fx - polynomial coefficients in bernstein basis
%
% gx - polynomial coefficients in bernstein basis
%
% ux - quotient polynomial coefficients
%
% vx - quotient polynomial coefficients
%
% alpha -
%
% theta -
%
% dx - initial approximation of GCD(fx,gx)
%
% t - Calculated degree of GCD
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Global Variables

global PLOT_GRAPHS
global MAX_ERROR_APF
global MAX_ITERATIONS_APF

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
m = size(fx,1) - 1;

% Get degree of polynomial g
n = size(gx,1) - 1;

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

% Get HCG
[HCG,H1C1G,H2C2G] = BuildHCG(uw,vw,m,n,t);

% Build HCG with respect to theta
[HCG_wrt_theta,H1C1G_wrt_theta,H2C2G_wrt_theta] = ...
    BuildHCG(Partial_uw_wrt_theta,Partial_vw_wrt_theta,m,n,t);

%Build the RHS vector b = [fw,alpha.*gw]
bk = [fw ; alpha(ite).*gw];

% Get initial residual
g = bk-((HCG)*dw);

% Set some initial values
residual(ite)   = norm(g);
p               = zeros(m+1,1);
q               = zeros(n+1,1);

% Values of LHS Perturbations
% Obtain structured perturbations sw of fw, and tw of gw
sw   = p.*(theta(ite).^vecm);
tw   = q.*(theta(ite).^vecn);


% Set initial values for the iterative process
z1x = zeros(length(uw),1);
z2x = zeros(length(vw),1);

% Construct the coefficient matrix in the equation that defines
% the constraint for the LSE problem.
HYk         = BuildHYQ(dx,m,n,theta(ite));


% % Build the matrix C given by Hz Hp Hq Halpha Htheta1 Htheta2 
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

C_temp      = ...
    [
    H_p,             zeros(m+1,n+1), zeros(m+1,1), H_theta1;...
    zeros(n+1,m+1),  H_q,            H_alpha,       H_theta2...
    ];

C       = [H_z , C_temp];

E       = eye(2*m+2*n-2*t+6);

fnew    = zeros(2*m+2*n-2*t+6,1);


ek = bk;

% Get the condition 
condition(ite) = norm(g)/norm(ek);

startpoint = [zk;p;q;alpha;theta];

yy = startpoint;

% Start the iterative procedure for the solution of the LSE problem.

while condition(ite) > (MAX_ERROR_APF) && ite < MAX_ITERATIONS_APF
    
    % Use the QR decomposition to solve the LSE problem.
    % min |y-p| subject to Cy=q
    y = LSE(E,fnew,C,g);
    
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
    [HCG,H1C1G,H2C2G] = BuildHCG(uw,vw,m,n,t);
     
    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G with respect to
    % theta
    [HCG_wrt_theta,H1C1G_wrt_theta, H2C2G_wrt_theta] = BuildHCG(uw_wrt_theta,vw_wrt_theta,m,n,t);
    
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
    [HEG,H1E1G,H2E2G] = BuildHCG(z1w,z2w,m,n,t);
    
    
    % Calculate Partial derivatives of Matrices H_{1}E_{1}(z1)G and
    % H_{2}E_{2}(z2)G with respect to theta
    [HEG_wrt_theta,H1E1G_wrt_theta,H2E2G_wrt_theta] = BuildHCG(z1w_wrt_theta,z2w_wrt_theta,m,n,t);
     
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
    
    % Build Matrix C
    
    % Calculate H_z
    HYk = BuildHYQ(dx,m,n,theta(ite));
    
    % Build H_z
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
    [HCEG,H1C1E1G,H2C2E2G] = BuildHCG(uw+z1w,vw+z2w,m,n,t);
    
    % Calculate the new residual
    g = [fw + sw ;(alpha(ite)*(gw + tw))]-((HCEG)*dw);
    
    % Calculate the new right hand vector.
    ek = [fw + sw ;(alpha(ite)*(gw + tw))];
    
    
    % update vector of residual
    residual(ite) = norm(g);
    
    % Update Condition scalar.
    condition(ite) = norm(g)/norm(ek);
    
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


if ite == MAX_ITERATIONS_APF
    fprintf('APF failed to converge after %i iterations\n',ite)
    PostAPF_dx = dx;
    PostAPF_ux = ux;
    PostAPF_vx = vx;
    PostAPF_theta = initial_theta;
    
    % Edit 20/07/2015
    PostAPF_fx = fx;
    PostAPF_gx = gx;
    figure()
    hold on
    plot(log10(condition))
    hold off
    return
    
end

try
    switch PLOT_GRAPHS
        case 'y'
            % Plot the normalised residuals res_ux, res_vx, res_uw and res_vw.
            plotgraphs3(res_ux,res_vx,res_uw,res_vw);
            
            % Write out the number of iterations required and plot the values of
            % alpha, theta and the residual.
            plotgraphs4(alpha,theta,residual);
        case 'n'
        otherwise
            error('err')
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



switch PLOT_GRAPHS
    case 'y'
        
        plot_residuals = 'y';
        plot_thetas = 'y';
        plot_betas = 'y';
        
        switch plot_residuals
            case 'y'
                figure('name','APF - Residuals')
                title('plotting residual')
                hold on
                plot(1:1:length(residual),(residual));
        end
        switch plot_thetas
            case 'y'
                figure('name','APF - thetas')
                title('plotting theta')
                hold on
                plot(1:1:length(theta),(theta));
        end
        switch plot_betas
            case 'y'
                figure('name','APF - alpha')
                title('plotting beta');
                hold on
                plot(1:1:length(alpha),(alpha));
        end
    case 'n'
    otherwise
        error('err')
end

end




