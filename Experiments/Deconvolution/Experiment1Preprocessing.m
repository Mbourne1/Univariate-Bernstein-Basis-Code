function [] = Experiment1Preprocessing(ex_num)

close all;
clc;


noise = 1e-4;

% Good Examples
% 8, 9, 10, 11,



bool_preproc = false;
o_Deconvolution(ex_num, noise, bool_preproc)


bool_preproc = true;
o_Deconvolution(ex_num, noise, bool_preproc)

end