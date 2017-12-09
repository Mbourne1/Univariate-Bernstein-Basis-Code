function [] = Experiment0(ex_num, bool_preproc)
%
%
%
% >> Experiment0('1')

% Examples

% 1 -
% 2 -
% 3 -
% 4 -
% 5 -
% 6 - Too small
%
%

%close all; clc;




emin = 1e-4;
emax = 1e-4;



switch bool_preproc
    case true
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
    case false
        mean_method = 'None';
        bool_alpha_theta = false;
end



apf_method = 'None';
sylvester_build_method = 'DTQ';
rank_revealing_metric = 'Minimum Singular Values';



%deconvolution_bool = false;


low_rank_approx_method = 'Standard SNTLN';
%low_rank_approx_method = 'None';

deconvolution_method_hx = 'Batch Constrained With STLN';
%deconvolution_method_hx = 'Batch Constrained';
%deconvolution_method_hx = 'Separate';


deconvolution_method_wx = 'Batch';
%deconvolution_method_wx = 'Separate'

deconvolution_preproc = true;



o_roots_Univariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, sylvester_build_method, ...
    rank_revealing_metric, deconvolution_method_hx, ...
    deconvolution_method_wx, deconvolution_preproc)

%SavePlots(ex_num)

%sameaxes()

end


function [] = SavePlots(ex_num)
%
% % Inputs
%
% ex_num : (String)
%
% str : (String)
%
% sylvester_build_method

directory_name = strcat('Univariate_Roots/Example',(ex_num),'/');

mkdir(directory_name)
h = get(0,'children');

%myFileName = strcat(sylvester_build_method, '_', str);
myFileName = 'Figure';
nPlots = length(h);

for i = 1 : 1 : nPlots
    
    %saveas(h(i), [directory_name num2str(length(h) + 1 - i)], 'fig');
    saveas(h(i), [directory_name strcat(myFileName, num2str(nPlots-i))], 'fig');
    saveas(h(i), [directory_name strcat(myFileName, num2str(nPlots-i))], 'eps');
    saveas(h(i), [directory_name strcat(myFileName, num2str(nPlots-i))], 'png');
    
end


end