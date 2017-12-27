% Get optimal mean method


close all; clc; 

ex_num = '8';

[fx_exact, gx_exact, dx_exact, ux_exact, vx_exact] = Examples_GCD(ex_num);

fx = fx_exact;
gx = gx_exact;

m = GetDegree(fx);
n = GetDegree(gx);
k = min(m,n) - 5;

AM_fx = GetArithmeticMean(fx, n - k);
AM_gx = GetArithmeticMean(gx, m - k);

GM_fx = GetGeometricMeanMatlabMethod(fx, n - k);
GM_gx = GetGeometricMeanMatlabMethod(gx, m - k);


fx_gm = fx ./ GM_fx;
gx_gm = gx ./ GM_gx;

fx_am = fx ./ AM_fx;
gx_am = gx ./ AM_gx;

% Build the kth subresultant

    S1 = BuildSubresultant_2Polys(fx, gx, k);
    S2 = BuildSubresultant_2Polys(fx_gm, gx_gm, k);
    S3 = BuildSubresultant_2Polys(fx_am, gx_am, k);
    
    bool_logs = true;
    if bool_logs == true
        S1 = real(log10(S1));
        S2 = real(log10(S2));
        S3 = real(log10(S3));
    end
    
    figure()
    subplot(3,1,1)
    h1 = heatmap(S1)
    
    subplot(3,1,2)
    h2 = heatmap(S2)
    
    subplot(3,1,3)
    h3 = heatmap(S3)