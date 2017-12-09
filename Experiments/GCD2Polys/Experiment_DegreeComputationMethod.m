function [] = Experiment_DegreeComputationMethod(ex_num, emax, bool_preproc)

% Variable : Preprocessing

% Good Examples

% Example 11 - 1e-10 , 1e-12


close all; clc;

%emin = 1e-10;
%emax = 1e-8;
%emin = emax * (10^(-2));
emin = emax;

switch bool_preproc
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
    case true
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
        
end

% Constants
low_rank_approx_method = 'None';
apf_method = 'None';
sylvester_build_method = 'DTQ';

%'R1 Row Diagonals'
%'R1 Row Norms'
%'Residuals'


arrRankRevealingMetric = {'Normalised Minimum Singular Values', 'Normalised R1 Row Diagonals', 'R1 Row Norms'};

for i = 1 : 1 : length(arrRankRevealingMetric)
    
    rank_revealing_metric = arrRankRevealingMetric{i};
    
    o_gcd_Univariate_2Polys(ex_num, emin, emax, mean_method, ...
        bool_alpha_theta, low_rank_approx_method, apf_method, ...
        sylvester_build_method, rank_revealing_metric) ;

end


end