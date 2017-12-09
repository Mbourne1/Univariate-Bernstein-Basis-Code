function [] = Experiment1SylvesterFormat_2Polys(ex_num)
% This example, we want to see the computation of the degree of the GCD
% with different forms of the sylvester subresultant matrix. 
%
% Variable : Sylvester Build Type
%
%
%
% >> Experiment1SylvesterFormat_2Polys('1')


% Bad Examples
%
% Without Preprocessing
%
% 20 - Too large without preprocessing
%



close all;
clc; 

% Constants 
el = 1e-14; 
eu = 1e-14;


mean_method = 'None';
bool_preproc = false;
 
%mean_method = 'Geometric Mean Matlab Method';
%bool_preproc = true;

rank_revealing_metric = 'Minimum Singular Values';

arrSylvesterFormat = {'T', 'DT','TQ', 'DTQ'};
low_rank_approx_method = 'None';


for i = 1:1: length(arrSylvesterFormat)

    sylvesterFormat = arrSylvesterFormat{i};
    
    
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