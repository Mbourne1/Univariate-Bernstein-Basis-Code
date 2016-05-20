function [fx_n,gx_n,alpha,theta] = LowRankApproximation(fx_n,gx_n,alpha,theta,t,opt_col,gm_fx,gm_gx)
%
%
% Inputs.
%
% fx_n : Normalized coefficients of polynomial f(x)
%
% gx_n : Normalized coefficients of polynomial g(x)
%
% alpha : Optimal value of alph
%
% theta : Optimal value of theta
%
% t : Degree of GCD of f(x) and g(x)
%
% opt_col : Index of optimal column for removal from S_{t}(f,g)
%
% gm_fx : Mean of entries of f(x) in T_{n-t}(f)
%
% gm_gx : Mean of entries of g(x) in T_{m-t}(g)


global SETTINGS

switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    case 'Standard STLN'
        
        fw = GetWithThetas(fx_n,theta);
        gw = GetWithThetas(gx_n,theta);
        
        % Performe structured Total Least Norm
        [fw,a_gw] = STLN(fw,alpha.*gw,t,opt_col);
        
        % Get f(x) and g(x) from low rank approximation.
        fx_n = GetWithoutThetas(fw,theta);
        
        gw = a_gw ./ alpha;
        gx_n = GetWithoutThetas(gw,theta);
        
        
    case 'Standard SNTLN' % Structured Non-Linear Total Least Norm
        
        % Perform Structured non-linear total least norm
        [fx_n,gx_n,alpha,theta,~] = SNTLN(fx_n,gx_n,alpha,theta,t,opt_col);
        
        
    case 'Root Specific SNTLN'
            
        [fx_n,gx_n,alpha,theta,~] = ...
            SNTLN_Roots(fx_n,gx_n,alpha,theta,t,opt_col,gm_fx,gm_gx);
        
    case 'None'
        
    otherwise
        error('SETTINGS.LOW_RANK_APPROXIMATION_METHOD must be valid')
end

end