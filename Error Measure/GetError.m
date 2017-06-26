
function [error] = GetError(name, fx_calc, fx_exact)
% Used in Computing the distance between the GCD and the computed GCD in
% the GCD finding method.
%
% Get distance between f(x) and the calulated f(x)
%
% Get the angle between the two vectors
% angle = dot(f_calc,f_exact) ./ (norm(f_calc) * norm(f_exact));
% angle_error = 1 - angle;
% fprintf('\tCalculated angle error : %8.2e \n', angle_error)

fx_calc_n  = NormaliseVector(fx_calc);
fx_calc_n = fx_calc_n ./ fx_calc_n(1);

fx_exact_n = NormaliseVector(fx_exact);
fx_exact_n = fx_exact_n ./ fx_exact_n(1);

% Calculate relative errors in outputs
absolute_error_fx = abs(fx_calc_n - fx_exact_n);

% Print Relative Error
relative_error_fx = norm( absolute_error_fx ./ fx_exact_n);
fprintf('\tCalculated relative error %s : %8.2e \n ',name, relative_error_fx);

% Print norm of absolute error
error = norm(absolute_error_fx);
fprintf('\tCalculated error %s : %8.2e \n', name, error);



end