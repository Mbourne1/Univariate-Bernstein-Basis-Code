% This experiment looks at the singular values of each sylvester
% subresultant format for 3 permutations of the same gcd problem


close all; clc;

emin = 1e-9;
emax = 1e-9;

% arrEx_num = {'11'}; emin = 1e-9; emax = 1e-9;

% Good Examples

% Bad Examples
% 9 - too small


mean_method = 'Geometric Mean Matlab Method';
bool_alpha_theta = true;

%mean_method = 'Geometric Mean Matlab Method';
%bool_alpha_theta = true;

low_rank_approx_method = 'None';
apf_method = 'None';

rank_revealing_metric = 'Minimum Singular Values';
bool_reorder_polys = false;


%arrSylvesterFormat = {'T All Equations', 'DT All Equations', 'TQ All Equations', 'DTQ All Equations'};
arrSylvesterFormat = {'T', 'DT', 'TQ', 'DTQ'};


ex_num = '13';




ex_num_variant = strcat(ex_num,'a');

for i1 = 1 : 1 : length(arrSylvesterFormat)
    sylvester_build_method = arrSylvesterFormat{i1};
    
    
    
    o_gcd_Univariate_3Polys(ex_num_variant, emin, emax, mean_method, ...
        bool_alpha_theta, low_rank_approx_method, apf_method,...
        sylvester_build_method, rank_revealing_metric, bool_reorder_polys)
    
    
    
end


