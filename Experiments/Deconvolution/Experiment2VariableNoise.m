function [] = Experiment2VariableNoise(ex_num)
% Perform a series of deconvolutions as found in the factorisation
% algorithm, for a variety of noise levels.
%
% % Inputs
%
% ex_num : (String) Example number


close all;
clc;


% Create an array of noise levels, where noise is added to the coefficients
% of the set of polynomials f_{i}(x)
arrNoise = {1e-6, 1e-8, 1e-10, 1e-12};

% Set whether to include or exclude preprocessing of the polynomials
% f_{i}(x)
bool_preproc = false;

% For each noise level, perform the deconvolution
for i = 1:1:length(arrNoise)
    
    noise = arrNoise{i};
    
    o_Deconvolution(ex_num, noise, bool_preproc)
    
    
end

end