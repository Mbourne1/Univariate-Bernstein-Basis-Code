function [] = Experiment(ex_num)

close all; clc;
% Get coefficients of f(x,y) and g(x,y)

[fx, gx] = Examples_GCD_FromCoefficients(ex_num);

% Get degree of polynomial f(x,y)
m = GetDegree(fx);

% Get degree of polynomial g(x,y)
n = GetDegree(gx);

k = 2;

global SETTINGS
SETTINGS.BOOL_LOG = false;
SETTINGS.SYLVESTER_BUILD_METHOD = 'DTQ';

% Get Geometric mean
gm_fx = GetGeometricMeanMatlabMethod(fx, n - k);
gm_gx = GetGeometricMeanMatlabMethod(gx, m - k);


% Get Arithmetic Mean
am_fx = GetArithmeticMean(fx, n - k);
am_gx = GetArithmeticMean(gx, m - k);


fx_n_gm = fx ./ gm_fx;
gx_n_gm = gx ./ gm_gx;


fx_n_am = fx ./ am_fx;
gx_n_am = gx ./ am_gx;


plot_logs = true;

if plot_logs
   fx_n_gm = log10( fx_n_gm);
   fx_n_am = log10( fx_n_am);
   fx = log10(fx);
   
   gx_n_gm = log10( gx_n_gm);
   gx_n_am = log10( gx_n_am);
   gx = log10(gx);
end


figure_name = 'fx and fx_n';
figure('name',figure_name)
hold on
x_vec = 0:1:m;
plot(x_vec, fx_n_gm, '-s', 'DisplayName','\bar{f}(x)_{gm}')
plot(x_vec, fx_n_am, '-s', 'DisplayName','\bar{f}(x)_{am}')
plot(x_vec, fx, '-s', 'DisplayName','f(x)')
legend(gca,'show');
hold off

figure_name = 'gx and gx_n';
figure('name', figure_name)
hold on
x_vec = 0:1:n;
plot(x_vec, gx_n_gm, '-s', 'DisplayName','\bar{f}(x)_{gm}')
plot(x_vec, gx_n_am, '-s', 'DisplayName','\bar{f}(x)_{am}')
plot(x_vec, gx, '-s', 'DisplayName','f(x)')
legend(gca,'show');
hold off



end



