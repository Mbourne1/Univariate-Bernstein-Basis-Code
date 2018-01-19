global SETTINGS

LineBreakLarge();
fprintf('PARAMETERS: \n')
fprintf('\t EXAMPLE NUMBER : %s \n', SETTINGS.EX_NUM)
fprintf('\t MINIMUM NOISE : %2.4e \n', emin)
fprintf('\t MAXIMUM NOISE : %2.4e \n', emax)
fprintf('\t MEAN METHOD : %s \n', SETTINGS.MEAN_METHOD)
fprintf('\t BOOL ALPHA THETA : %s \n', num2str(SETTINGS.BOOL_ALPHA_THETA))
fprintf('\t LOW RANK APPROXIMATION METHOD : %s \n', SETTINGS.LOW_RANK_APPROXIMATION_METHOD);
fprintf('\t APF METHOD : %s \n ', SETTINGS.APF_METHOD)
fprintf('\t SYLVESTER MATRIX VARIANT : %s \n', SETTINGS.SYLVESTER_MATRIX_VARIANT)
fprintf('\t RANK REVEALING METIC : %s \n', SETTINGS.RANK_REVEALING_METRIC)
fprintf('\t DECONVOLUTION METHOD HX : %s \n', SETTINGS.DECONVOLUTION_METHOD_HX)
fprintf('\t DECONVOLUTION METHOD WX : %s \n', SETTINGS.DECONVOLUTION_METHOD_WX)
fprintf('\t BOOL LOG: %s \n', num2str(SETTINGS.BOOL_LOG))

LineBreakLarge();