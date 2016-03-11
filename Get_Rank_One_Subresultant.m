function [t] = Get_Rank_One_Subresultant()
% Given that only one subresultant exists (S_{1}(f,g) and min(m,n) = 1,
% the two polynomials f(x) and g(x) are either coprime, or have a gcd of
% degree 1, the GCD must be a multiple of either f(x) or g(x).

global THRESHOLD


% If the minimum degree is one, only one subresultant exists.

fprintf('Only One Sylvester Subresultant Exists \n')
fprintf('min(m,n) = 1 \n')
fprintf('Degree of GCD is either one or zero \n')

[~,R1] = qr(Sk);

%% Analysis of Diagonal entries of R1

% Get the diagonals of R1
vDiags_R1 = diags(R1);

% Get the vector of changes in diagonal values of R1.
vDelta_Diags_R1 = abs(diff(log10(vDiags_R1)));

% Get the maximum change in the diagonal entries of R1
maxDelta_diags_R1 = max(vDelta_Diags_R1);

%% Analysis of Row sums of R1

% Get the vector of row sums of R1
vRowSum_R1 = sum(R1,2);

% Get the changes in row sums
vDelta_RowSum_R1 = abs(diff(log10(vRowSum_R1)));

% Get the maximum change ni row sum
max_delta_log_rowsum_R1 = max(log10(vDelta_RowSum_R1));

%%
if max_delta_log_rowsum_R1 > THRESHOLD
    % The maximum change in row sum (delta) is significant
    % Subresultant is rank deficient
    % Set degree of GCD = 1.
    
    t = 1;
    
    return;
    
else
    % The maximum change in row sum (delta) is insignificant
    % Subresultant is of full rank
    % Set degree of GCD = 0
    
    % Set Degree of GCD to zero
    t = 0;
    return;
end