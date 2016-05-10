function [fx_n,gx_n,alpha,theta] = LowRankApproximation(fx_n,gx_n,alpha,theta,t,opt_col)


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