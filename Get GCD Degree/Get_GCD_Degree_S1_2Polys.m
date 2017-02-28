function [t, alpha, theta, GM_fx, GM_gx] = Get_GCD_Degree_One_S1_2Polys(fx, gx)
% Compute the degree of the GCD of two polynomials by computing the
% singular values of S(f,g).
%
% % Inputs
%
% [fx, gx] : Coefficients of polynomials f(x) and g(x)
%
% % Outputs
%
% t
%
% alpha
%
% theta
%
% GM_fx
%
% GM_gx


k = 1;
i = 1;

% Get degree of the polynomial f(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Preprocess polynomials f(x) and g(x) as they appear in S_{1}(f,g)
[GM_fx, GM_gx, alpha, theta] = Preprocess_2Polys(fx,gx,k);

% Divide f(x) and g(x) by geometric means
fx_n = fx./GM_fx(i);
gx_n = gx./GM_gx(i);

% Construct the kth subresultant matrix S_{k}(f(\theta),g(\theta))
fw = GetWithThetas(fx_n, theta(i));
gw = GetWithThetas(gx_n, theta(i));

% Build the k-th subresultant
Sk = BuildSubresultant_2Polys(fw, alpha(i).*gw, k);

vSingularValues = svd(Sk);

global SETTINGS
if (SETTINGS.PLOT_GRAPHS)
    figure_name = sprintf([mfilename sprintf('Singular Values of %s',SETTINGS.SYLVESTER_BUILD_METHOD)]);
    figure('name',figure_name)
    hold on
    plot(log10(vSingularValues),'-s');
    hold off
end

vMetric = vSingularValues;

% Get the change in the ratios from one subresultant to the next.
vDeltaMetric = abs(diff(log10(vMetric)));

% Get the maximum change in rowsum ratio and its index
[maxDelta, indexMaxDelta] = max(vDeltaMetric);


t = (m+n) - indexMaxDelta;

end
