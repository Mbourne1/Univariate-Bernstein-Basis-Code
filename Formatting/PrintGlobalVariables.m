global SETTINGS
LineBreakLarge();
fprintf('PARAMETERS: \n')
fprintf('\t Min noise : %2.4e \n', emin)
fprintf('\t Max noise : %2.4e \n', emax)
fprintf('\t Mean Method : %s \n', SETTINGS.MEAN_METHOD)
fprintf('\t Alpha Theta Boolean : %s \n', num2str(SETTINGS.BOOL_ALPHA_THETA))
fprintf('\t Low Rank Approximation Method : %s \n', SETTINGS.LOW_RANK_APPROXIMATION_METHOD);
fprintf('\t APF Method : %s \n ', SETTINGS.APF_METHOD)
fprintf('\t Log Bool: %s \n', num2str(SETTINGS.BOOL_LOG))
fprintf('\t Sylvester Build Method : %s \n', SETTINGS.SYLVESTER_BUILD_METHOD)
LineBreakLarge();