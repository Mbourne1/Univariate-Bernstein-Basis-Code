fprintf('PARAMETERS:\n\n')
fprintf('\tExample Number: %i \n',ex_num)
fprintf('\tmin noise : %2.4e \n\tmax noise : %2.4e \n',emin,emax)
fprintf('INPUT VARIABLES')
fprintf('\n\tSNTLN : %s \n',...
    bool_sntln);
fprintf('\tAPF : %s \n ',BOOL_APF)
fprintf('\tDENOM : %s \n',BOOL_DENOM_SYL)
fprintf('\tPREPROC : %s \n',BOOL_PREPROC)
fprintf('\tLOG: %s \n',BOOL_LOG)
fprintf('\tQ : %s \n',BOOL_Q)
fprintf('\tSNTLN Derivative Constraint: %s \n',BOOL_SNTLN_ROOTS)
fprintf('')
fprintf('--------------------------------------------------------------------------- \n')