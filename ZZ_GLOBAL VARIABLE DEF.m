% BOOL_APF - (Boolean)
%    1 :- Apply Structured Perturbations to the Approximate Polynomial
%    Factorisation, before obtaining GCD d_{k}.
%    0 :- Don't Include
%    Structured Perturbations.

% BOOL_Q - (Boolean) Consists of binomial coefficients of coprime
% polynomials.
%   1 :- Include the binomial coefficients from the null space in the
%    Sylvester Matrix.
%   0 :- Exclude the binomial coefficients from the null space in the
%    Sylvester Matrix.

% BOOL_SNTLN - (Boolean)
%    1 :- Include Structured Perturbations in Sylvester Matrix S_{k} before
%    calculating quotients u_{k}, v_{k} and GCD d_{k}.
%    0 :- Don't Include Structured Perturbations.

% BOOL_PREPROC - (Boolean)
% It has been shown that the inclusion of preprocessors Geometric mean,
% scaling by alpha, change of independent variable, yield improved results.
%   1 :- Include Preprocessors.
%   0 :- Exclude Preprocessors.