function [] = CoefficientPlotting_2Polys
% Plot the coefficients of the polynomials of a set of examples and save
% the figures.


% Create array of example numbers
ex_num_arr = {'18','19','20','21','22'};

% For each example, plot and save the coefficients of f(x) and g(x)
for i = 1 : 1 : length(ex_num_arr)
    
    ex_num = ex_num_arr{i};
    plotCoefficients(ex_num);
    
end



end


function[] = plotCoefficients(ex_num)
% Plot coefficients of a given example
%
% % Inputs
%
% ex_num : (String) Example number


[fx_exact, gx_exact, ~, ~, ~] = Examples_GCD(ex_num);


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


