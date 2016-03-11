function C = BuildC(u,v,t)
% Build the Matric C
%
% Inputs.
%
% u : Polynomial u(x) or u(w)
%
% v : Polynomial v(x) or v(w)
%
% t : Degree of the GCD.

% Build partition of C corresponding to polynomial u(x)
C1 = BuildC1(u,t);

% Build partition of C corresponding to polynomial v(x)
C2 = BuildC1(v,t);

% Form C
C = [C1 ; C2];
end