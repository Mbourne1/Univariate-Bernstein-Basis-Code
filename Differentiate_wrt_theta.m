function Partial_fw_wrt_theta = Differentiate_wrt_theta(fw,theta)
% Differentiate f(\omega,\theta) with respect to theta
%
% % Inputs
%
% fw : Polynomial f(\omega)
%
% theta : Optimal value of \theta


bool_diff_method = 'power basis';
bool_diff_method = 'Bernstein basis';

switch bool_diff_method
    case 'power basis'
        
        Partial_fw_wrt_theta = Differentiate_wrt_theta_powerbasis(fw,theta);
        
    case 'Bernstein basis'
        
        Partial_fw_wrt_theta = Differentiate_wrt_theta_Bernsteinbasis(fw,theta);
        
    otherwise
        
        error('err')
        
end


end