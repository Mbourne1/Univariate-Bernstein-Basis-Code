function [residual,idx_col] = GetMinDistance(Sk)
% Given a Sylvester matrix remove each column in turn, and compute the
% distance between it and a vector which lies in the column space of the
% remaining columns.
% 
% Inputs.
%
% Sk : Sylvester subresultant matrix S_{k}(f,g)
%
% Outputs.
%
% residual : Residual obtained by removing the column ck. r = ck - Ax
%
% idx_col : Index of column removed where residual is minimal.

% Get the number of columns in the Sylvester matrix S_{k}(f,g)
[nCols] = size(Sk,2);

% Initialise a vector to store residuals
residual_vec = zeros(nCols,1);

for i = 1:1:nCols
    
    % Get Ak, Sk with ck removed
    Ak = Sk;
    Ak(:,i) = [];
    
    % Get the column ck
    ck = Sk(:,i);
    
    % Get the solution x
    x_ls = SolveAx_b(Ak,ck);
    
    % Get the residual
    residual_vec(i) = norm(ck - (Ak*x_ls));
    
end

% Get the minimal residual and the index of the corresponding column ck
[residual,idx_col] = min(residual_vec);

end