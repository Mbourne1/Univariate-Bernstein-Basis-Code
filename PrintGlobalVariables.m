global LOW_RANK_APPROXIMATION_METHOD
global APF_METHOD
global BOOL_DENOM_SYL
global BOOL_ALPHA_THETA
global MEAN_METHOD

global BOOL_LOG
global BOOL_Q

fprintf('PARAMETERS:\n\n')
fprintf('\tmin noise : %2.4e \n\tmax noise : %2.4e \n',emin,emax)
fprintf('INPUT VARIABLES\n')

fprintf('\t DENOM : %s \n',BOOL_DENOM_SYL)
fprintf('\t MEAN METHOD : %s \n',MEAN_METHOD)
fprintf('\t ALPHA_THETA : %s \n',BOOL_ALPHA_THETA)
fprintf('\t Low Rank Approximation Method : %s \n',LOW_RANK_APPROXIMATION_METHOD);
fprintf('\t APF Method : %s \n ',APF_METHOD)
fprintf('\t LOG: %s \n',BOOL_LOG)
fprintf('\t Q : %s \n',BOOL_Q)
fprintf('')
fprintf('--------------------------------------------------------------------------- \n')