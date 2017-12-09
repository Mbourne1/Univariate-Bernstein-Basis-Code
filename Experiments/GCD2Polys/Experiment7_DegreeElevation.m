function [] = Experiment7_DegreeElevation(ex_num, bool_preproc, p, q)


%ex_num = '1';
emin = 1e-12;
emax = 1e-12;

switch bool_preproc
    case true
        
        mean_method = 'Geometric Mean Matlab Method';
        bool_alpha_theta = true;
    case false
        
        mean_method = 'None';
        bool_alpha_theta = false;
end

low_rank_approx_method = 'None';
apf_method = 'None';
Sylvester_Build_Method = 'DTQ';
rank_revealing_metric = 'Minimum Singular Values';
%p = 5;
%q = 5;

o_gcd_Univariate_2Polys_DegreeElevationTest(ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, ...
    Sylvester_Build_Method, rank_revealing_metric, p, q)

end