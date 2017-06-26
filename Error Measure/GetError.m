
function [error] = GetError(name, fx_calc,fx_exact)
% Used in Computing the distance between the GCD and the computed GCD in
% the GCD finding method.
%
% Get distance between f(x) and the calulated f(x)
%
% Get the angle between the two vectors
% angle = dot(f_calc,f_exact) ./ (norm(f_calc) * norm(f_exact));
% angle_error = 1 - angle;
% fprintf('\tCalculated angle error : %8.2e \n', angle_error)

fx_calc  = NormaliseVector(fx_calc);

fx_exact = NormaliseVector(fx_exact);

% Calculate relative errors in outputs
rel_error_fx = norm(abs(fx_calc - fx_exact) ./ fx_exact);

fprintf('\tCalculated relative error %s : %8.2e \n ',name, rel_error_fx);

error = norm(abs(fx_calc - fx_exact) );

fprintf('\tCalculated error %s : %8.2e \n', name,error);



end