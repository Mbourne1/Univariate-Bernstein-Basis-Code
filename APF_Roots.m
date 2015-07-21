function [PostAPF_fx, PostAPF_gx, PostAPF_dx, PostAPF_uk, PostAPF_vk, PostAPF_theta] = ...
    APF_Roots(fx,ux,vx,initial_theta,dx,t)
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

% Set initial value of theta.
theta(ite) = initial_theta;

% Set the initial value of alpha.
%alpha(ite) = initial_alpha;

% Get degree of polynomial f
m = length(fx) - 1;

% Initialise some useful vectors
vecm  = (0:1:m)';
vecmk = (0:1:m-t)';
veck  = (0:1:t)';

% Convert f and g to modified bernstein basis, excluding binomial
% coefficient
fw_n = fx.*(theta(ite).^vecm);


% Convert u and v to modified bernstein basis, excluding binomial
% coefficient
uw = ux.*(theta(ite).^vecmk);

% Build the C(u);C(v) matrix the old method
H1C1G = BuildHCGPart(uw,t);
HCG = H1C1G;

% Convert d to modified bernstein basis, excluding binomial coefficient
dw = dx.*(theta(ite).^veck);


% Initialise S, such that sk = S * p_{t}, where p_{t} is a vector of
% perturbations applied to u
S = (diag(theta(ite).^vecm));

% Initialise zk - Structured perturbations of u and v
zk = zeros(m-t+1,1);

% Get partial derivatives of fw and gw with respect to theta
Partial_fw_wrt_theta = vecm.*(fw_n ./ theta(ite));

% Get partial derivative of dw with respect to theta
Partial_dw_wrt_theta = veck.*(dw./theta(ite));

% Get partial derivatives of uw and vw with respect to theta.
Partial_uw_wrt_theta = vecmk .*(uw./theta(ite));

% Build H_{t}C(f,g)G_{t}
switch Bool_APFBuildMethod
    case 0
        H1 = BuildH1(m);
        C1 = BuildC1Partition(uw,t);
        G  = BuildG(t);
        H1C1G = H1*C1*G;
    case 1
        
        % Edit 02/06/2015 17:00:00
        % Change the way HCG matrix is built
        H1C1G       = BuildHCGPart(uw,t);
        
end



% Build HCG with respect to theta
switch Bool_APFBuildMethod
    case 0 % Build Naive method, calculate H C and G independently
        C1_wrt_theta    = BuildC1Partition(Partial_uw_wrt_theta,t);
        H1C1G_wrt_theta = H1*C1_wrt_theta*G;

    case 1 % Build the matrix HCG
        H1C1G_wrt_theta     = BuildHCGPart(Partial_uw_wrt_theta,t);   
end

%Build the RHS vector b = [fw,alpha.*gw]
bk = fw_n;

% Get initial residual
rk = bk-((HCG)*dw);

% Set some initial values
gnew            = rk;               % The initial residual vector
residual(ite)   = norm(gnew);       % the initial normalized residual
p               = zeros(m+1,1);     % Perturbations to polynomial f

% Values of LHS Perturbations

% Set initial values for the iterative process
z1x = zeros(length(uw),1);


% Construct the coefficient matrix in the equation that defines
% the constraint for the LSE problem.

HYk        = BuildHYQ_Roots(dx,m,theta(ite));

uw_binom = zeros(length(uw),1);
for i = 0:1:length(uw)-1
    uw_binom(i+1) = uw(i+1) .* nchoosek(m-t,i);
end


H_z         = HYk;

H_p         = (-1)*S;

H_theta1    = -Partial_fw_wrt_theta+...
    (H1C1G_wrt_theta * dw)+...
    (H1C1G * Partial_dw_wrt_theta);

C_temp      = [ H_p,            H_theta1];

C       = [H_z , C_temp];

E       = eye(2*m-t+3);

fnew    = zeros(2*m-t+3,1);

ek = bk;

condition = norm(rk)/norm(ek);

startpoint = [zk;p;theta];

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
    
    % % obtain the small changes.

    delta_zk    = y(1:m-t+1);
    delta_pk    = y(m-t+2:2*m-t+2);
    delta_theta = y(end);
    
    % % Update variables zk, pk, qk, beta, theta
    
    zk          = zk + delta_zk;
    p           = p  + delta_pk;
    
    theta(ite)  = theta(ite-1) + delta_theta;
    
    % Update the iterative value of f and g
    fw = fx.*(theta(ite).^vecm);

    % Update matrices S and T
    S = diag(theta(ite).^vecm);

    % Update vectors of thetas
    ok1 = theta(ite).^vecmk;
    ok3 = theta(ite).^veck;
    
    % Update the iterative value of GCD dw
    dw = dx.*(ok3);
    
    % Update the iterative value of quotients uw and vw
    uw = ux.*(ok1);

    % Obtain partial derivatives of uw and vw with respect to theta
    uw_wrt_theta = vecmk .* (uw ./ theta(ite));

    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G
    switch Bool_APFBuildMethod
        case 0
            C1      = BuildC1Partition(uw,t);
            H1C1G   = H1*C1*G;
        case 1
            H1C1G   = BuildHCGPart(uw,t); 
    end

    % Build Matrices H_{1}C_{1}(u)G and H_{2}C_{2}(v)G with respect to
    % theta
    switch Bool_APFBuildMethod
        case 0
            C1_wrt_theta    = BuildC1Partition(uw_wrt_theta,t);
            H1C1G_wrt_theta = H1*C1_wrt_theta*G;

        case 1   
            H1C1G_wrt_theta = BuildHCGPart(uw_wrt_theta,t);
      
    end
    
    % Get perturbation vector, and seperate in to perturbations of f,
    % z1 and perturbations of g, z2
    z1x = zk(1:m-t+1);

    % Obtain z1 and z2 in the modified bernstein basis.
    z1w = z1x .* ok1;

    % obtain partial derivatives of z1 and z2 with respect to theta
    z1w_wrt_theta = vecmk .* (z1w ./ theta(ite));
    
    % Build Matrices H_{1}E_{1}(z1)G and H_{2}E_{2}(z2)G
    switch Bool_APFBuildMethod
        case 0
            E1 = BuildC1Partition(z1w,t);
            H1E1G = H1*E1*G;

        case 1 
            H1E1G = BuildHCGPart(z1w,t);

    end
    
    % Calculate Partial derivatives of Matrices H_{1}E_{1}(z1)G and
    % H_{2}E_{2}(z2)G with respect to theta
    
    switch Bool_APFBuildMethod
        case 0
            E1_wrt_theta = BuildC1Partition(z1w_wrt_theta,t);
            H1E1G_wrt_theta = H1 * E1_wrt_theta * G;
 
        case 1
            H1E1G_wrt_theta = BuildHCGPart(z1w_wrt_theta,t);
            
    end
    
    
    % Obtain structured perturbations sw of fw, and tw of gw
    sw   = p.*(theta(ite).^vecm);
   
    % Calculate partial derivatives of sw and tw with respect to theta
    
    Partial_s_wrt_theta     = vecm.*(p.*(theta(ite).^(vecm-1)));
    
    % Calculate partial derivatives of fw and gw with respect to theta
    Partial_fw_wrt_theta    = vecm.*(fx.*(theta(ite).^(vecm-1)));
   
    % Calculate partial derivative of dw with respect to theta
    Partial_dw_wrt_theta    = veck.*(dx.*(theta(ite).^(veck-1)));
    
    % Calculate H_z
       
    HYk_new         = BuildHYQ_Roots(dx,m,theta(ite));
    
    H_z = HYk_new;
    
    % Calculate H_p
    H_p = (-1)*S;
    
    % Calculate H_theta1
    H_theta1 = -(Partial_fw_wrt_theta + Partial_s_wrt_theta)+...
        ((H1C1G_wrt_theta)*dw)+...
        ((H1E1G_wrt_theta)*dw)+...
        ((H1C1G)*Partial_dw_wrt_theta)+...
        ((H1E1G)*Partial_dw_wrt_theta);

    C_temp = [  H_p, H_theta1];
    
    % Build Matrix C
    C = [H_z , C_temp];
    
    % Calculate Matrix H(C+E)G
    switch Bool_APFBuildMethod
        case 0
            H1C1E1G = H1 * (C1+E1) * G;
            HCEG = H1C1E1G;
        case 1
            H1C1E1G = BuildHCGPart(uw+z1w,t);
            HCEG    = H1C1E1G;
    end
    
    % Calculate the new residual
    rk = fw + sw - ((HCEG)*dw);
    
    % Calculate the new right hand vector.
    ek = fw + sw;
    
    % Update gnew - used in LSE Problem
    gnew = rk;
    
    % update vector of residual
    residual(ite) = norm(gnew);
    
    % Update Condition scalar.
    condition = norm(rk)/norm(ek);
    
    % Update fnew
    fnew = -(yy - startpoint);
    
    
end

if ite == max_iterations 
    fprintf('APF failed to converge after %i iterations\n',ite)
    PostAPF_dx = dx; 
    PostAPF_uk = ux;
    PostAPF_vk = vx;
    PostAPF_theta = initial_theta;
    
    % Edit 20/07/2015
    PostAPF_fx = fx;
    PostAPF_gx = gx;
    return 
end


% Display number of iterations
fprintf('--------------------------------------------------------------------------- \n')
fprintf('Iterations over approximate polynomial factorisation : %i \n', ite+1);
fprintf('--------------------------------------------------------------------------- \n')

% Update values of quotients u and v,
PostAPF_uk = ux + z1x;
PostAPF_vk = vx;

% Update value of theta
PostAPF_theta = theta(ite);

% Update value of common divisor dx
PostAPF_dx = dw./(theta(ite).^veck);

% update value of the input polynomial f
PostAPF_fx = (fw + sw) ./ (theta(ite) .^ vecm);

% update values of the input polynomial g, which is now calculated by
% derivative constraint.
PostAPF_gx = zeros(m,1);
% for i = 0:1:m-1
%    PostAPF_gx(i+1) = m.*(PostAPF_fx(i+2) - PostAPF_fx(i+1)); 
% end




bool_plotgraphs = 1;

switch bool_plotgraphs
    case 1
        
        plot_residuals = 1;
        plot_thetas = 1;
        
        switch plot_residuals
            case 1
                figure(100)
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

end

end









