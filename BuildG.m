function G = BuildG(t)
% Build G - A Diagonal matrix which corresponds to the binomial
% coefficients of the GCD t.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs.


% t :   Degree of GCD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize a zero vector
G = zeros(t+1,1);

% For each coefficient of the GCD, get is corresponding binomial 
% coefficient.
for i = 0:1:t
    G(i+1) = nchoosek(t,i);
end

% Diagonalise G
G = diag(G);

end
