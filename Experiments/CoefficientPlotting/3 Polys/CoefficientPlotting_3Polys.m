function [] = CoefficientPlotting_3Polys


ex_num_arr = {'1', '2','3','4','5','6','7','8','9','10','11','12','14','15','16'};



for i = 1 : 1 : length(ex_num_arr)

    ex_num = ex_num_arr{i};
    variant = 'a';
    
    ex_num = strcat(ex_num,variant);
    
    plotCoefficients(ex_num);
    
end



end


function[] = plotCoefficients(ex_num)

[fx_exact, gx_exact, hx_exact, dx_exact, ux_exact, vx_exact, wx_exact] = ...
    Examples_GCD_3Polys(ex_num);


PlotCoefficients({fx_exact,gx_exact, hx_exact},...
    {'$f(x)$', '$g(x)$', '$h(x)$'},...
    {'-s', '-o','-*'} ...
    )

file_name = strcat('Example_',ex_num);
formattype = 'epsc';
saveas(gcf,file_name,formattype)

formattype = 'png';
saveas(gcf,file_name,formattype)

formattype = 'jpeg';
saveas(gcf,file_name,formattype)
end


