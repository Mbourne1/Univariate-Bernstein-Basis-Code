function H1 = BuildH1(m)
% Build the matrix H1, used in BuildH().
% 
% Inputs
%
% m : Degree of polynomial.
%

H = diag(1./GetBinomials(m));

end