function Partial_fw_wrt_theta = Differentiate_wrt_theta(fw,th)
% Differentiate f(\omega,\theta) with respect to theta

% Get the degree of polynomial f(\omega)
m = GetDegree(fw);

% Initialise the multiplication vector.
vecm = (0:1:m)';

% Get the partial derivative of f with respect to theta
Partial_fw_wrt_theta = vecm .*fw ./ th;

end