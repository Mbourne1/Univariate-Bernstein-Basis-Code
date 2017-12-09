

close all;
clc;

% Good Examples
% 8, 9, 10, 11, 

ex_num = '12';
arrNoise = {1e-6, 1e-8, 1e-10, 1e-12};
bool_preproc = false;

for i = 1:1:length(arrNoise)

    noise = arrNoise{i};
    o_Deconvolution(ex_num, noise, bool_preproc)


end