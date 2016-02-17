function [roots_Bb] = o_roots_multroot(fx)
% Given a vector of coefficients of a polynomail in Bernstein form. Return
% the set of roots given by zheng MultRoots() function.
%
% % Inputs.
%
%   fx : Column vector of coefficients of the polynomial f(x) in Bernstein
%   form.
%
% % Outputs.
%
%   roots_Bb : Roots of f(x) as calculated by the zheng MultRoots function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

addpath 'Root Finding Methods/multroot/multroot'

% Get the degree of polynomial f(x) 
[r,~] = size(fx);
m = r-1;

% Build the vector of corresponding binomial coefficients
Bi_m = GetBinomials(m);

roots_calc = multroot(fx.*Bi_m);
    

% convert roots wrt t, Bernstein basis
roots_Bb = [1- roots_calc(:,1)./(1+roots_calc(:,1)) roots_calc(:,2)];

% Printout roots to screen
PrintoutRoots('MULTROOTS METHOD', roots_Bb);




end

