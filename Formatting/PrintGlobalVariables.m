global SETTINGS

fprintf('PARAMETERS:\n\n')
fprintf('\tmin noise : %2.4e \n\tmax noise : %2.4e \n', emin, emax)
fprintf('INPUT VARIABLES\n')


fprintf('\t MEAN METHOD : %s \n', SETTINGS.MEAN_METHOD)
fprintf('\t ALPHA_THETA : %s \n', num2str(SETTINGS.BOOL_ALPHA_THETA))
fprintf('\t Low Rank Approximation Method : %s \n', SETTINGS.LOW_RANK_APPROXIMATION_METHOD);
fprintf('\t APF Method : %s \n ', SETTINGS.APF_METHOD)
fprintf('\t LOG: %s \n', SETTINGS.BOOL_LOG)
fprintf('\t Sylvester Build Method : %s \n', SETTINGS.SYLVESTER_BUILD_METHOD)
fprintf('')
fprintf('--------------------------------------------------------------------------- \n')