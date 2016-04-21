global LOW_RANK_APPROXIMATION_METHOD
global APF_METHOD
global BOOL_DENOM_SYL
global BOOL_ALPHA_THETA
global MEAN_METHOD

global BOOL_LOG
global BOOL_Q

fprintf('PARAMETERS:\n\n')
fprintf('\tmin noise : %2.4e \n\tmax noise : %2.4e \n',emin,emax)
fprintf('INPUT VARIABLES')
fprintf('\tDENOM : %s \n',BOOL_DENOM_SYL)
fprintf('\t MEAN : %s \n',MEAN_METHOD)
fprintf('\tALPHA_THETA : %s \n',BOOL_ALPHA_THETA)
fprintf('\n\tLow Rank Approximation Method : %s \n',...
    LOW_RANK_APPROXIMATION_METHOD);
fprintf('\tAPF Method : %s \n ',APF_METHOD)

fprintf('\tLOG: %s \n',BOOL_LOG)
fprintf('\tQ : %s \n',BOOL_Q)
fprintf('')
fprintf('--------------------------------------------------------------------------- \n')