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






% Variables

arrSylvesterFormat = {'DTQ'};

%arr_ex_num_variant = {strcat(ex_num,'a'), strcat(ex_num,'b'), strcat(ex_num,'c')} ;
arr_ex_num_variant = {strcat(ex_num,'c')};

arr_noise = {1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4};






% Set number of equations to define the sylvester matrix
% 2 : Sylvester matrix has a 2 x 3 partitioned structure
% 3 : Sylvester matrix has a 3 x 3 partitioned structure

nEquations = '2';

% For each subresultant matrix variant
for i1 = 1 : 1 : length(arrSylvesterFormat)
    
    % For each example number variant 'a', 'b' or 'c'
    for i2 = 1 : 1 : length(arr_ex_num_variant)
        
        
        for i3 = 1 : 1 : length(arr_noise)
            
            emin = arr_noise{i3};
            emax = arr_noise{i3};
            
            
            % Set subresultant matrix variant
            sylvester_matrix_variant = arrSylvesterFormat{i1};
            
            % Set example number
            ex_num_variant = arr_ex_num_variant{i2};
            
            o_gcd_Univariate_3Polys(ex_num_variant, emin, emax, mean_method, ...
                bool_alpha_theta, low_rank_approx_method, apf_method,...
                sylvester_matrix_variant, nEquations, rank_revealing_metric)
            
            %SavePlots(ex_num, 'fgh', sylvester_matrix_variant)
            %close all; clc;
        end
    end
    
end








end

function [] = SavePlots(ex_num, str, sylvester_matrix_variant)
%
% % Inputs
%
% ex_num : (String)
%
% str : (String)
%
% sylvester_matrix_variant : (String)

directory_name = strcat('Example',(ex_num),'/Figures_',str,'/');

mkdir(directory_name)
h = get(0,'children');

myFileName = strcat(sylvester_matrix_variant, '_', str);

for i=2:1:2
    
    %saveas(h(i), [directory_name num2str(length(h) + 1 - i)], 'fig');
    saveas(h(i), [directory_name myFileName], 'fig');
    saveas(h(i), [directory_name myFileName], 'eps');
    saveas(h(i), [directory_name myFileName], 'png');
    
end


end





