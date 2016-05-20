function [PostAPF_fx, PostAPF_gx, PostAPF_dx, PostAPF_ux, PostAPF_vx, PostAPF_theta] = ...
    APF(fx,gx,ux,vx,i_alpha,i_theta,dx,t)
% Refine Approximate Polynomial Factorisation by APF
%
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
% i_alpha - input alpha
%
% i_theta - input theta
%
% dx - initial approximation of GCD(fx,gx)
%
% t - Calculated degree of GCD
%
%

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
veck    = (0:1:t)';

% Convert f and g to modified bernstein basis, excluding binomial
% coefficient
fw = GetWithThetas(fx,th);
gw = GetWithThetas(gx,th);

% Convert u and v to modified bernstein basis, excluding binomial
% coefficient
uw = GetWithThetas(ux,th);
vw = GetWithThetas(vx,th);

% Convert d to modified bernstein basis, excluding binomial coefficient
dw = GetWithThetas(dx,th);

% Initialise S, such that sk = S * pt
S = (diag(th(ite).^vecm));

% Initialise T - Matrix such that tk = T * qt
T = (diag(th(ite).^vecn));

% Initialise zk - Structured perturbations of u and v
zk = zeros(m+n-2*t+2,1);

% Get partial derivatives of fw and gw with respect to theta
fw_wrt_theta = Differentiate_wrt_theta(fw,th(ite));
gw_wrt_theta = Differentiate_wrt_theta(gw,th(ite));

% Get partial derivative of dw with respect to theta
dw_wrt_theta = Differentiate_wrt_theta(dw,th(ite));

% Get partial derivatives of uw and vw with respect to theta.
Partial_uw_wrt_theta = Differentiate_wrt_theta(uw,th(ite));
Partial_vw_wrt_theta = Differentiate_wrt_theta(vw,th(ite));

% Get HCG
[HCG,H1C1G,H2C2G] = BuildHCG(uw,vw,m,n,t);

% Build HCG with respect to theta
[~,H1C1G_wrt_theta,H2C2G_wrt_theta] = ...
    BuildHCG(Partial_uw_wrt_theta,Partial_vw_wrt_theta,m,n,t);

%Build the RHS vector b = [fw,alpha.*gw]
bk = [fw ; alpha(ite).*gw];

% Get initial residual
res_vec = bk-((HCG)*dw);

% Set some initial values
residual(ite)   = norm(res_vec);
p               = zeros(m+1,1);
q               = zeros(n+1,1);

% Values of LHS Perturbations
% Obtain structured perturbations sw of fw, and tw of gw
sw = GetWithThetas(p,th);
tw = GetWithThetas(q,th);

% Set initial values for the iterative process
z1x = zeros(length(uw),1);
z2x = zeros(length(vw),1);

% Construct the coefficient matrix in the equation that defines
% the constraint for the LSE problem.
HYk         = BuildHYQ(dx,m,n,th(ite));


% % Build the matrix C given by Hz Hp Hq Halpha Htheta1 Htheta2
H_z         = HYk;

H_p         = (-1)*S;

H_q         = -(alpha(ite))*T;

H_alpha     = -(gw + tw);

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

E       = eye(2*m+2*n-2*t+6);

fnew    = zeros(2*m+2*n-2*t+6,1);


ek = bk;

% Get the condition
condition(ite) = norm(res_vec)/norm(ek);

startpoint = [zk;p;q;alpha;th];

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
    th(ite)     = th(ite-1) + delta_theta;
    
    % Update the iterative value of f and g
    fw = GetWithThetas(fx,th(ite));
    gw = GetWithThetas(gx,th(ite));
    
    % Update matrices S and T
    S = diag(th(ite).^vecm);
    T = diag(th(ite).^vecn);
    
    % Update the iterative value of GCD dw
    dw = GetWithThetas(dw,th(ite));
    
    % Update the iterative value of quotients uw and vw
    uw = GetWithThetas(ux,th(ite));
    vw = GetWithThetas(vx,th(ite));
    
    % Obtain partial derivatives of uw and vw with respect to theta
    uw_wrt_theta = Differentiate_wrt_theta(uw,th(ite));
    vw_wrt_theta = Differentiate_wrt_theta(vw,th(ite));
    
    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G
    [~,H1C1G,H2C2G] = BuildHCG(uw,vw,m,n,t);
    
    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G with respect to
    % theta
    [~,H1C1G_wrt_theta, H2C2G_wrt_theta] = BuildHCG(uw_wrt_theta,vw_wrt_theta,m,n,t);
    
    % Get perturbation vector, and seperate in to perturbations of f,
    % z1 and perturbations of g, z2
    z1x = zk(1:m-t+1);
    z2x = zk(m-t+2:m+n-2*t+2);
    
    % Obtain z1 and z2 in the modified bernstein basis.
    z1w = GetWithThetas(z1x,th(ite));
    z2w = GetWithThetas(z2x,th(ite));
    
    % obtain partial derivatives of z1 and z2 with respect to theta
    z2w_wrt_theta = Differentiate_wrt_theta(z2w,th(ite));
    z1w_wrt_theta = Differentiate_wrt_theta(z1w,th(ite));
    
    % Build Matrices H_{1}E_{1}(z1)G and H_{2}E_{2}(z2)G
    [~,H1E1G,H2E2G] = BuildHCG(z1w,z2w,m,n,t);
    
    % Calculate Partial derivatives of Matrices H_{1}E_{1}(z1)G and
    % H_{2}E_{2}(z2)G with respect to theta
    [~,H1E1G_wrt_theta,H2E2G_wrt_theta] = BuildHCG(z1w_wrt_theta,z2w_wrt_theta,m,n,t);
    
    % Obtain structured perturbations sw of fw, and tw of gw
    sw = GetWithThetas(p,th(ite));
    tw = GetWithThetas(q,th(ite));
    
    % Calculate partial derivatives of sw and tw with respect to theta
    s_wrt_theta = Differentiate_wrt_theta(sw,th(ite));
    t_wrt_theta = Differentiate_wrt_theta(tw,th(ite));
    
    % Calculate partial derivatives of fw and gw with respect to theta
    fw_wrt_theta = Differentiate_wrt_theta(fw,th(ite));
    gw_wrt_theta = Differentiate_wrt_theta(gx,th(ite));
    
    % Calculate partial derivative of dw with respect to theta
    dw_wrt_theta = Differentiate_wrt_theta(dw,th(ite));
    
    % Build Matrix C
    
    % Calculate H_z
    HYk = BuildHYQ(dx,m,n,th(ite));
    
    % Build H_z
    H_z = HYk;
    
    % Calculate H_p
    H_p = (-1)*S;
    
    % Calculate H_q
    H_q = -(alpha(ite))*T;
    
    % Calculate H_beta
    H_alpha = -(gw + tw);
    
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
    [HCEG,~,~] = BuildHCG(uw+z1w,vw+z2w,m,n,t);
    
    % Calculate the new residual
    res_vec = [fw + sw ;(alpha(ite)*(gw + tw))]-((HCEG)*dw);
    
    % Calculate the new right hand vector.
    ek = [fw + sw ;(alpha(ite)*(gw + tw))];
    
    
    % update vector of residual
    residual(ite) = norm(res_vec);
    
    % Update Condition scalar.
    condition(ite) = norm(res_vec)/norm(ek);
    
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
    
    fx_p = GetWithoutThetas(fw,th(ite));
    sx_p = GetWithoutThetas(sw,th(ite));
    gx_p = GetWithoutThetas(gw,th(ite));
    tx_p = GetWithoutThetas(tw,th(ite));
    ukx = ux + z1x;
    vkx = vx + z2x;
    dkx = dw./(th(ite).^veck);
    
    [res_ux(ite),res_vx(ite)] = Term_Criterion_APF(fx_p,gx_p,sx_p,...
        tx_p,ukx,vkx,dkx,t,1.0);
    
    
end


if ite == SETTINGS.MAX_ITERATIONS_APF
    fprintf('APF failed to converge after %i iterations\n',ite)
    PostAPF_dx = dx;
    PostAPF_ux = ux;
    PostAPF_vx = vx;
    PostAPF_theta = i_theta;
    
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
catch
    fprintf('Can not perform plotgraphs3 and plotgraphs4, not enough iterations performed\n')
end
% END EDIT 01/06/2015 16:34:00


% Display number of iterations
fprintf('Iterations over approximate polynomial factorisation : %i \n', ite+1);

% Update values of quotients u and v,
PostAPF_ux = ux + z1x;
PostAPF_vx = vx + z2x;

% Update value of theta
PostAPF_theta = th(ite);

% Update value of common divisor dx
PostAPF_dx = GetWithoutThetas(dw,th(ite));

% Edit 20/07/2015
PostAPF_fx = GetWithoutThetas((fw + sw),th(ite));
PostAPF_gx = GetWithoutThetas((gw + tw),th(ite));



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




