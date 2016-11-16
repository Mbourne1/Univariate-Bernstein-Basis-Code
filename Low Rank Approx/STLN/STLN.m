function [fx_lr,gx_lr,ux_lr,vx_lr] = STLN(fx,gx,k,idx_col)
% Perform STLN with no refinement of alpha or theta to compute the low rank
% approximation of the Sylvester subresultant matrix S_{k}(f,g)
%
% Inputs.
%
% fx : Coefficients of polynomial f(x)
%
% gx : coefficients of polynomial g(x)
%
% k : Degree of GCD d(x)
%
% idx_col : Index of optimal column to be removed from S_{t}(f,g)
%
% Outputs.
%
% fx_lr : Coefficients of f(x) after addition of structured perturbations
%
% gx_lr : Coefficients of g(x) after addition of strucutred perturbations
%
% ux_lr : Coefficients of u(x)
%
% vx_lr : Coefficients of v(x)

global SETTINGS

% Get the degree of polynomial f(x)
m = GetDegree(fx);

% Get the degree of polynomial g(x)
n = GetDegree(gx);

% Initialise the vector of perturbations zf(x)
z_fx = zeros(m+1,1);

% Initialise the vector of perturbations zg(x)
z_gx = zeros(n+1,1);

% Initialise the vector of perturbations z.
z = [z_fx ; z_gx];

% Build the t'th subresultant S_{t}(f,g)
DTQ = BuildDTQ(fx,gx,k);

% Build the matrix E_{t}(z)
DBQ = BuildDTQ(z_fx,z_gx,k);

% Get the index of the optimal colummn for removal
%[~,colIndex] = GetMinDistance(St);

% Get A_{t} the LHS matrix, equivalent to S_{t} with the optimal column
% removed.
Ak = DTQ;
Ak(:,idx_col) = [];

% Get c_{t} the removed column of S_{t} to form A_{t}.
ck = DTQ(:,idx_col);

% Get E_{t}, the matrix of strucured perturbations corresponding to A_{t}.
Bk = DBQ;
Bk(:,idx_col) = [];

% Get h_{t}, the vector of strucutred perturbations corresponding to c_{t}
hk = DBQ(:,idx_col);

% Build DP
DPQ = BuildDP_SNTLN(m,n,1,1,idx_col,k);

% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
xk = SolveAx_b(Ak+Bk,ck+hk);

% Get residual vector
res_vec = (ck + hk) - (Ak+Bk)*xk;

% Get the vector x with a zero included in the x_ls solution.
x = [xk(1:idx_col-1) ; 0 ; xk(idx_col:end)];

% Build the matrix Y_{t}
DYQ = BuildDYQ_SNTLN(x,m,n,k,1,1);

% Build the Matrx C for LSE problem

H_z = DYQ - DPQ;
H_x = Ak + Bk;

C = [H_z H_x];

% Build the identity matrix E.
E = blkdiag(eye(m+n+2),zeros(m+n-2*k+1,m+n-2*k+1));


% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    z;...
    xk;
    ];


% Initialise yy the vector of accumulated perturbations.
yy = start_point;

% Set the initial value of vector p to be zero
f = -(yy-start_point);


% Initialise the iteration counter
ite = 1;

% Set the termination criterion.
condition(ite) = norm(res_vec)./norm(ck);


while condition(ite) >  SETTINGS.MAX_ERROR_SNTLN &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
    
    % increment iteration number
    ite = ite + 1;
    
    % Get small changes in vector y
    y_lse = LSE(E,f,C,res_vec);
    
    % add small changes to cummulative changes
    yy = yy + y_lse;
    
    % obtain the small changes in z and x
    delta_zk        = y_lse(1:m+n+2,1);
    delta_xk        = y_lse((m+n+3):(2*m+2*n-2*k+3),1);
    
    % Update z and x
    z = z + delta_zk;
    xk = xk + delta_xk;
    
    % Split z into z_f and z_g
    z_fx = z(1:m+1);
    z_gx = z(m+2:end);
    
    % Build the matrix E = D^{-1} * B(zf,zg) * Q
    DBQ = BuildDTQ(z_fx,z_gx,k);
    
    % Build the matrix Bt = DBQ with opt column removed.
    Bk = DBQ;
    Bk(:,idx_col) = [];
    
    % Get h_{t}, the optimal column removed from BDQ
    hk = DBQ(:,idx_col);
    
    % Get the vector x
    x = [xk(1:idx_col-1) ; 0 ; xk(idx_col:end)];
    
    % Seperate the component parts of x into x_v and x_u, where x_v is an
    % approximation of v(x) and x_u is an approximation u(x).
    DYQ = BuildDYQ_SNTLN(x,m,n,k,1,1);
    
    % Get updated residual vector
    res_vec = (ck + hk) - ((Ak + Bk)*xk);
    
    % Update the matrix C
    H_z = DYQ - DPQ;
    H_x = Ak + Bk;
    C = [H_z H_x];
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    % Update the termination criterion.
    condition(ite) = norm(res_vec) ./ norm(ck + hk) ;
    
end


% %
% Get polynomials for output.

% Get polynomial f(x) + \delta f(x)
fx_lr = fx + z_fx;

% Get polynomial g(x) + \delta g(x)
gx_lr = gx + z_gx;

% %
% Get u(x) and v(x) from x_ls
x = [xk(1:idx_col-1) ; -1 ; xk(idx_col:end)];

% Get polynomial u(x)
nCoeff_vx = n-k+1;
vx_lr = x(1:nCoeff_vx);

% Get polynomial v(x)
ux_lr = -1.* x(nCoeff_vx + 1:end) ;

% % 
% Plot Graphs
Plot_STLN()

LineBreakLarge()
fprintf([mfilename ' : ' sprintf('Required number of iterations : %i \n',ite)]);
LineBreakLarge()
end


