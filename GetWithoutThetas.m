function fw = GetWithoutThetas(fx,th)
% Take the coefficients of polynomial f(\omega,\theta), divide each 
% coefficient a_{i} by \theta^{i} to obtain f(x)

m = GetDegree(fx);

% Multiply coefficients of f(x) a_{i} by \theta^{i}

fw = fx ./ (th.^(0:1:m)');


end