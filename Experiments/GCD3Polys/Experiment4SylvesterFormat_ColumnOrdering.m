function [] = Experiment4SylvesterFormat_ColumnOrdering(ex_num, bool_preproc)
%
% This experiment looks at the singular values of each sylvester
% subresultant format for 3 permutations of the same gcd problem
%
% >> Experiment4SylvesterFormat_ColumnOrdering('7')


% % Good Examples
%
% ex_num = 7, emin = 1e-7; emax = 1e-3;


close all; clc;

emin = 1e-9;
emax = 1e-9;

% arrEx_num = {'11'}; emin = 1e-9; emax = 1e-9;

% Good Examples

switch bool_preproc
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
        
    otherwise
        error('err')
        
end
low_rank_approx_method = 'None';
apf_method = 'None';

rank_revealing_metric = 'Minimum Singular Values';


arrSylvesterFormat = {'DTQ'};

arr_ex_num_variant = {strcat(ex_num,'a'), strcat(ex_num,'b'), strcat(ex_num,'c')} ;





% Set number of equations to define the sylvester matrix
% 2 : Sylvester matrix has a 2 x 3 partitioned structure
% 3 : Sylvester matrix has a 3 x 3 partitioned structure

nEquations = '2';

for i1 = 1 : 1 : length(arrSylvesterFormat)
    
    for i2 = 1 : 1 : length(arr_ex_num_variant)
        
        sylvester_build_method = arrSylvesterFormat{i1};
        
        ex_num_variant = arr_ex_num_variant{i2};
        
        o_gcd_Univariate_3Polys(ex_num_variant, emin, emax, mean_method, ...
            bool_alpha_theta, low_rank_approx_method, apf_method,...
            sylvester_build_method, nEquations, rank_revealing_metric)
        
        %SavePlots(ex_num, 'fgh', sylvester_build_method)
        %close all; clc;
    end
    
end

nEquations = '3';
ex_num_variant = strcat(ex_num,'a');

o_gcd_Univariate_3Polys(ex_num_variant, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method,...
    sylvester_build_method, nEquations, rank_revealing_metric)







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

directory_name = strcat('Example',(ex_num),'/Figures_',str,'/');

mkdir(directory_name)
h = get(0,'children');

myFileName = strcat(sylvester_build_method, '_', str);

for i=2:1:2
    
    %saveas(h(i), [directory_name num2str(length(h) + 1 - i)], 'fig');
    saveas(h(i), [directory_name myFileName], 'fig');
    saveas(h(i), [directory_name myFileName], 'eps');
    saveas(h(i), [directory_name myFileName], 'png');
    
end


end





