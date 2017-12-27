function [] = Experiment1SylvesterFormat_2Polys(ex_num)
% This example, we want to see the computation of the degree of the GCD
% with different forms of the sylvester subresultant matrix.
%
% % Inputs
%
% ex_num : (String) Example number
%
%
% Variable : Sylvester Build Type
%
%
% >> Experiment1SylvesterFormat_2Polys('1')

close all;
clc;

% Constants
el = 1e-14;
eu = 1e-14;


% Determine whether preprocessing is included
bool_preprocessing = true;
if (bool_preprocessing == true)
    mean_method = 'Geometric Mean Matlab Method';
    bool_preproc = true;
else
    mean_method = 'None';
    bool_preproc = false; 
end


% Constants
rank_revealing_metric = 'Minimum Singular Values';
arrSylvesterVariants = {'T', 'DT','TQ', 'DTQ'};
low_rank_approx_method = 'None';


for i = 1 : 1: length(arrSylvesterVariants)
    
    sylvesterFormat = arrSylvesterVariants{i};
    
    
    o_gcd_Univariate_2Polys(ex_num, el, eu, mean_method, bool_preproc, ...
        low_rank_approx_method, 'None', sylvesterFormat, rank_revealing_metric);
    
    %SavePlots(ex_num, '',sylvesterFormat);
    %close all; clc;
end


end



function [] = SavePlots(ex_num, str, sylvester_build_method)
%
% % Inputs
%
% ex_num : (String)
%
% str : (String)
%
% sylvester_build_method

directory_name = strcat('UnivariateSylvesterFormatFigures/Example',(ex_num),'/Figures_',str,'/');

mkdir(directory_name)
h = get(0,'children');

myFileName = strcat(sylvester_build_method, '_', str);

for i=2:1:2
    try
        %saveas(h(i), [directory_name num2str(length(h) + 1 - i)], 'fig');
        saveas(h(i), [directory_name myFileName], 'fig');
        saveas(h(i), [directory_name myFileName], 'eps');
        saveas(h(i), [directory_name myFileName], 'png');
    catch
        error('err')
    end
end


end