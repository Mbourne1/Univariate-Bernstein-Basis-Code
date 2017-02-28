


% Set variables
noise_level = 1e-8;
ex_num_arr = {'1','2','3','4','5','6','7','8'};
ex_num_arr = {'9','10'};
sylvester_method_arr = {'DTQ','TQ','DT','T'};


% Global variables
global SETTINGS
SETTINGS.MEAN_METHOD = 'Geometric Mean Matlab Method';
SETTINGS.BOOL_ALPHA_THETA = 'y';
SETTINGS.SEED = 1024;
SETTINGS.GCD_COEFFICIENT_METHOD = 'ux and vx';
SETTINGS.APF_BUILD_METHOD = 'Standard';



err_ux_unproc = zeros(length(ex_num_arr),length(sylvester_method_arr));
err_vx_unproc = zeros(length(ex_num_arr),length(sylvester_method_arr));
err_dx_unproc = zeros(length(ex_num_arr),length(sylvester_method_arr));
err_ux_preproc = zeros(length(ex_num_arr),length(sylvester_method_arr));
err_vx_preproc = zeros(length(ex_num_arr),length(sylvester_method_arr));
err_dx_preproc = zeros(length(ex_num_arr),length(sylvester_method_arr));
opt_col_mat = zeros(length(ex_num_arr),length(sylvester_method_arr));
opt_col_preproc_mat = zeros(length(ex_num_arr),length(sylvester_method_arr));
condition_mat = zeros(length(ex_num_arr),length(sylvester_method_arr));
condition_preproc_mat = zeros(length(ex_num_arr),length(sylvester_method_arr));

warning('OFF')
for i1 = 1:1:length(ex_num_arr)
    
    ex_num = ex_num_arr{i1};
    [fx_exact, gx_exact, dx_exact, ux_exact, vx_exact] = Examples_GCD(ex_num);
    fx = AddNoiseToPoly(fx_exact,noise_level);
    gx = AddNoiseToPoly(gx_exact,noise_level);
    
    
    
    
    m = GetDegree(fx_exact);
    n = GetDegree(gx_exact);
    t = GetDegree(dx_exact);
    
    
    for i2 = 1:1:length(sylvester_method_arr)
        
        SETTINGS.SYLVESTER_BUILD_METHOD = sylvester_method_arr{i2};
        Sk_unproc = BuildSubresultant_2Polys(fx,gx,t);
        
        [~,idx] = GetMinDistance(Sk_unproc);
        opt_col_mat(i1,i2) = idx;
        
        [Ak,ck] = RemoveSubresultantColumn(Sk_unproc,idx);
        condition_mat(i1,i2) = cond(Ak);
        
        [ux_unproc, vx_unproc] = GetQuotients_2Polys(fx,gx,t);
        err_ux_unproc(i1,i2) = GetError(ux_unproc,ux_exact);
        err_vx_unproc(i1,i2) = GetError(vx_unproc,vx_exact);
        
        dx_unproc = GetGCDCoefficients_2Polys(ux_unproc,vx_unproc,fx,gx,t,1,1);
        err_dx_unproc(i1,i2) = GetError(dx_unproc,dx_exact);
        
        
        % % % PREPROCESS
        
        [GM_fx, GM_gx, alpha, theta] = Preprocess_2Polys(fx,gx,t);
        fx_n = fx ./ GM_fx;
        gx_n = gx ./ GM_gx;
        fw = GetWithThetas(fx_n, theta);
        gw = GetWithThetas(gx_n, theta);
        
        
        Sk_preproc = BuildSubresultant_2Polys(fw,alpha.*gw,t);
        
        [~,idx_preproc] = GetMinDistance(Sk_preproc);
        opt_col_preproc_mat(i1,i2) = idx_preproc;
                      
        [Ak_preproc,ck_preproc] = RemoveSubresultantColumn(Sk_preproc,idx_preproc);
        condition_preproc_mat(i1,i2) = cond(Ak_preproc);
        
        [uw, vw] = GetQuotients_2Polys(fw, alpha.*gw, t);
        
        ux_preproc = GetWithoutThetas(uw,theta);
        vx_preproc = GetWithoutThetas(vw,theta);
        
        err_ux_preproc(i1,i2) = GetError(ux_preproc,ux_exact);
        err_vx_preproc(i1,i2) = GetError(vx_preproc,vx_exact);
        
        dx_preproc = GetGCDCoefficients_2Polys(ux_preproc, vx_preproc, fx_n , gx_n, t, alpha, theta);
        
        
        err_dx_preproc(i1,i2) = GetError(dx_preproc,dx_exact);
        
        
        %[dx_preproc./dx_preproc(1) dx_exact./dx_exact(1)]
        
    end

    
    
end

warning('ON')

for i = 1:1:length(ex_num_arr)
    LineBreakLarge()
    
    fprintf(sprintf('Example %s \n',ex_num_arr{i}))
    fprintf(sprintf('m = %i, n = %i, t = %i \n',m,n,t));
    fprintf('\n');
    
    fprintf('Index of Removed Column \n')
    disp([opt_col_mat(i,:); opt_col_preproc_mat(i,:)]);
    
    fprintf('Condition \n')
    disp([condition_mat(i,:) ; condition_preproc_mat(i,:)]);
    
    fprintf('Error u(x)\n')
    disp([err_ux_unproc(i,:) ; err_ux_preproc(i,:)]);
    
    fprintf('Error v(x)\n')
    disp([err_vx_unproc(i,:) ; err_vx_preproc(i,:)]);
    
    fprintf('Error d(x)\n')
    disp([err_dx_unproc(i,:) ; err_dx_preproc(i,:)]);
    
end

function err = GetError(fx,fx_exact)


error_measure = 'angle';
switch error_measure
    
    case 'angle'
        
        %error is given by cos \theta
        
        fx = fx./fx(1);
        fx_exact = fx_exact./fx_exact(1);
        
        err = dot(fx,fx_exact) ./ (norm(fx) * norm(fx_exact));
        
    case 'relative'
        fx = fx./fx(1);
        fx_exact = fx_exact./fx_exact(1);
        
        rel_difference = abs(fx - fx_exact) ./ fx_exact;
        err = norm(rel_difference);
        
    otherwise
        error('err')
end
end



