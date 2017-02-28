function [] = Untitled



% Test 1 - Change noise level and save all RRM plots

ex_num = '2';
emax_arr = {1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15,1e-16,1e-17};
mean_method = 'None';
bool_alpha_theta = 'n';
low_rank_approx_method = 'None';
apf_method = 'None';
Sylvester_Build_Method = 'DTQ';

for i = 1:1:length(emax_arr)
    emax = emax_arr{i};
    o_gcd_2Polys(ex_num,emax,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method,Sylvester_Build_Method)   
    
   
    
    
end

sameaxes()
SaveAllFigures()


% Test 2 - Change Build method and save all RRM plots.



% Constants



%o_gcd_2Polys(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method,Sylvester_Build_Method)






end