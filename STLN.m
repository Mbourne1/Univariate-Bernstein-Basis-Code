function [fx,gx] = STLN(fx,gx,t,colIndex)
% Perform STLN with no preprocessors

global MAX_ERROR_SNTLN
global MAX_ITERATIONS_SNTLN
global PLOT_GRAPHS

% Get degree of polynomial f(x)
m = size(fx,1) - 1;

% Get the degree of polynomial g(x)
n = size(gx,1) - 1;

% Initialise the vector of perturbations zf(x)
zf = zeros(m+1,1);

% Initialise the vector of perturbations zg(x)
zg = zeros(n+1,1);

% Initialise the vector of perturbations z.
z = [zf ; zg];

% Build the t'th subresultant S_{t}(f,g)
DTQ = BuildDTQ(fx,gx,1,t);

% Build the matrix E_{t}(z)
DBQ = BuildDTQ(zf,zg,1,t);

% Get the index of the optimal colummn for removal
%[~,colIndex] = GetMinDistance(St);

% Get A_{t} the LHS matrix, equivalent to S_{t} with the optimal column
% removed.
At = DTQ;
At(:,colIndex) = [];

% Get c_{t} the removed column of S_{t} to form A_{t}.
ct = DTQ(:,colIndex);

% Get E_{t}, the matrix of strucured perturbations corresponding to A_{t}.
Bt = DBQ;
Bt(:,colIndex) = [];

% Get h_{t}, the vector of strucutred perturbations corresponding to c_{t}
ht = DBQ(:,colIndex);

% Build Pt
Pt = BuildPt(colIndex,m,n,t);

% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
x_ls = SolveAx_b(At+Bt,ct+ht);

% Get residual vector
res_vec = (ct + ht) - (At+Bt)*x_ls;

% Get the vector x with a zero included in the x_ls solution.
x = [x_ls(1:colIndex-1) ; 0 ; x_ls(colIndex:end)];

% Build the matrix Y_{t}
Yt = BuildYt(x,m,n,t);

% Build the Matrx C for LSE problem

H_z = Yt - Pt;
H_x = At + Bt;

C = [H_z H_x];

% Build the identity matrix E.
E = eye(2*m+2*n-2*t+3);
%E = blkdiag(eye(m+n+2),zeros(m+n-2*t+1,m+n-2*t+1))


% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    z;...
    x_ls;
    ];

% Set the initial value of vector p to be zero
f = zeros(2*m+2*n-2*t+3,1);

% Initialise yy the vector of accumulated perturbations.
yy = start_point;

% Initialise the iteration counter
ite = 1;

% Set the termination criterion.
condition(ite) = norm(res_vec)./norm(ct);


while condition(ite) >  MAX_ERROR_SNTLN &&  ite < MAX_ITERATIONS_SNTLN
    
    % increment iteration number
    ite = ite + 1;
    
    % Get small changes in vector y
    y_lse = LSE(E,f,C,res_vec);
    
    % add small changes to cummulative changes
    yy = yy + y_lse;
    
    % obtain the small changes in z and x
    delta_zk        = y_lse(1:m+n+2,1);
    delta_xk        = y_lse((m+n+3):(2*m+2*n-2*t+3),1);
    
    % Update z and x
    z = z + delta_zk;
    x_ls = x_ls + delta_xk;
    
    % Split z into z_f and z_g
    zf = z(1:m+1);
    zg = z(m+2:end);
    
    % Build the matrix E = DBQ
    DBQ = BuildDTQ(zf,zg,1,t);
    
    % Build the matrix Bt = DBQ with opt column removed.
    Bt = DBQ;
    Bt(:,colIndex) = [];
    
    % Get h_{t}, the optimal column removed from BDQ
    ht = DBQ(:,colIndex);
    
    % Get the vector x_ls
    x = [x_ls(1:colIndex-1) ; 0 ; x_ls(colIndex:end)];
    
    % Build matrix Y_{t}
    Yt = BuildYt(x,m,n,t);
    
    % Get updated residual vector
    res_vec = (ct+ht) - ((At+Bt)*x_ls);
    
    % Update the matrix C
    H_z = Yt - Pt;
    H_x = At + Bt;
    C = [H_z H_x];
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    % Update the termination criterion.
    condition(ite) = norm(res_vec) ./ norm(ct+ht) ;
    
end

% update f(x)
fx = fx + zf;
gx = gx + zg;

switch PLOT_GRAPHS
    case 'y'
        figure('name','STLN - Residuals')
        hold on
        plot(log10(condition),'-s');
        hold off
    case 'n'
    otherwise
        error('err')
end

fprintf('Required number of iterations : %i \n',ite)

end


function Pt = BuildPt(idx_Col,m,n,t)
% Build the matrix P_{t}, where h_{t} = P_{t}z

if idx_Col <= n-t+1
    % First Partition
    j = idx_Col;
    
    % Get binomials corresponding to f(x)
    bi_m = GetBinomials(m);
    
    bi_denom = zeros(m+1,1);
    for k = 0:1:m
        bi_denom(k+1) = nchoosek(m+n-t,(j-1)+k);
    end
    
    
    G = bi_m .* nchoosek(n-t,j-1) ./ bi_denom;
    
    Pt = ...
        [
        zeros(j-1,m+1)      zeros(j-1,n+1);
        diag(G)        zeros(m+1,n+1);
        zeros(n-t-j+1,m+1)  zeros(n-t-j+1,n+1);
        ];
else
    % Second Partition
    j = idx_Col - (n-t+1);
    
    
    bi_n = GetBinomials(n);
    
    bi_denom = zeros(n+1,1);
    for k = 0:1:n
        bi_denom(k+1) = nchoosek(m+n-t,j+k-1);
    end
    
    G = (bi_n .* nchoosek(m-t,j-1)) ./ bi_denom;
    
    Pt = ...
        [
        zeros(j-1,m+1)      zeros(j-1,n+1);
        zeros(n+1,m+1)      diag(G) ;
        zeros(m-t-j+1,m+1)  zeros(m-t-j+1,n+1);
        ];
    
end

end

function Yt = BuildYt(x,m,n,t)
% Build the matrix Y_{t}, where D*Y_{t}*K*z = D*E_{t}*Q*x

xa = x(1:n-t+1);
xb = x(n-t+2:end);

% Get x with binomials
bi_m = GetBinomials(m);
bi_n = GetBinomials(n);

D = BuildD(m,n,t);

Y1 = BuildC1(xa,m);
Y2 = BuildC1(xb,n);

Yt = D*[Y1 Y2]*blkdiag(diag(bi_m),diag(bi_n));


end