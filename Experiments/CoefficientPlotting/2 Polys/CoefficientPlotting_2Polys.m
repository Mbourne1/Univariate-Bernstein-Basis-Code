function [] = CoefficientPlotting_2Polys


%ex_num_arr = {'1', '2','3','4','5','6','7','8','9','10'...
    %,'11','12','13','14','15','16','17'};

    ex_num_arr = {'18','19','20','21','22'};
    
for i = 1 : 1 : length(ex_num_arr)

    ex_num = ex_num_arr{i};
    
    plotCoefficients(ex_num);
    
end



end


function[] = plotCoefficients(ex_num)

[fx_exact, gx_exact, dx_exact, ux_exact, vx_exact] = Examples_GCD(ex_num);


PlotCoefficients({fx_exact,gx_exact},...
    {'$f(x)$', '$g(x)$'},...
    {'-s', '-o'} ...
    )

file_name = strcat('Example_',ex_num);
formattype = 'epsc';
saveas(gcf,file_name,formattype)

formattype = 'png';
saveas(gcf,file_name,formattype)

formattype = 'jpeg';
saveas(gcf,file_name,formattype)
end


