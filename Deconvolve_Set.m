function hi = Deconvolve_Set(set_g)
% Performs a series of d deconvolutions over a set of polynomials,
% where each polynomial g_{i} appears in two deconvolutions.
%
%
% % Inputs
%
% set_g :   set of input polynomials g(y) to be deconvolved. Each g_{i} has a
%           different number of elements, so set_g is a cell array.
%
% % Outputs
%
% h_{i} = g_{i-1}/g_{i}
%
%



% Set Deconvolution Method
%   single := Use standard non batch method
%   batch  := Use batch deconvolution
global SETTINGS

switch SETTINGS.DECONVOLVE_METHOD
    case 'Separate'
        % Deconvolve independent method
        hi = Deconvolve_Independent(set_g);
    case 'Batch'
        % Deconvolve Batch Method
        hi = Deconvolve_Batch(set_g);
    otherwise
        error('SETTINGS.DECONVOLVE_METHOD must be either (Separate) or (Batch)')
end




end
