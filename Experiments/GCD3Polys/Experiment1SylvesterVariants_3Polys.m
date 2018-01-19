function [] = Experiment1SylvesterVariants_3Polys(arrEx_num, bool_preproc)
%
% % Inputs
%
% arrEx_num : (Array of Strings) Array of example numbers
%
% % Outputs


% Set noise level

emin = 1e-8;
emax = 1e-8;


low_rank_approx_method = 'None';
apf_method = 'None';


rank_revealing_metric = 'Minimum Singular Values';


%arrSylvesterFormat = {'T', 'DT', 'TQ', 'DTQ'};
arrSylvesterFormat = {'DTQ'};

nEquations = '2';


switch bool_preproc
    
    case 1
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
    case 0
        mean_method = 'None';
        bool_alpha_theta = false;
end

for j = 1 : 1 : length(arrEx_num)
    
    
    ex_num = arrEx_num{j};
    
    
    for i1 = 1 : 1 : length(arrSylvesterFormat)
        
        sylvester_matrix_variant = arrSylvesterFormat{i1};

        %try
        o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, ...
            bool_alpha_theta, low_rank_approx_method, apf_method,...
            sylvester_matrix_variant, nEquations, rank_revealing_metric)
        %catch
        
        
        %end
        
        
    end
end

nEquations = '3';

o_gcd_Univariate_3Polys(ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method,...
    sylvester_matrix_variant, nEquations, rank_revealing_metric)
