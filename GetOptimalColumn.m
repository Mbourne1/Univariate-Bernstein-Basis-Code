function opt_col = GetOptimalColumn(Sk)

% Return the optimal column of S_{k}, such that its removal to form a 
% right hand vector ck and the remaining columns form matrix A_{k}, 
% where ck - Ak*x gives minimal residual.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                               Input:

% Sk :- Sylvester matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                              Output

% opt_col :- The optimal column c_{k}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get number of columns in Sylvester Matrix
    [~,cols] = size(Sk);

%Initialise vector to store residuals
    residuals_QR = zeros(1,cols);
    [~,n] = size(Sk);
    n = n-1;

% QR Decomposition of the Sylvester Matrix S_{k}
    [Qk,Rk] = qr(Sk);

% For each column i, remove column from subresultant S_{k} and perform QR
% decomposition on S_{k,i}, Alternatively use QRdelete on Sk, for each
% column index i.
    for i = 1:1:cols
        % Get the ith column of S_{k}
        ck = Sk(:,i);
        % Perform QR delete on the Sylvester Matrix such that QR now
        % represents the QR Decomposition of S_{k,i} where S_{k,i} is the
        % Sylvester Matrix with ith column removed.
        [Q,~] = qrdelete(Qk,Rk,i);
        cd = Q'*ck;
        d = cd(n+1:end,:);
        
        residuals_QR(i) = norm(d);

    end

%Obtain the column for which the residual is minimal.
    [~,opt_col] = min(log10(residuals_QR));

%fprintf('Optimal column for removal is %i\n\n',opt_col);


end