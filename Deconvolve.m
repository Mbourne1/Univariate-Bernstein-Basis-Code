function hi = Deconvolve(set_g)
% Performs a series of d deconvolutions over a set of polynomials,
% where each polynomial g_{i} appears in two deconvolutions.
%
%
% Inputs
%
%
% set_g :   set of input polynomials g(y) to be deconvolved. Each g_{i} has a
%           different number of elements, so set_g is a cell array.
%
% Outputs
%
%
% h_{i} = g_{i-1}/g_{i}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Set Deconvolution Method
%   single := Use standard non batch method
%   batch  := Use batch deconvolution
global BOOL_DECONVOLVE

switch BOOL_DECONVOLVE
    case 'Single'
        % Deconvolve independent method
        hi = Deconvolve_Independent(set_g);
    case 'Batch'
        % Deconvolve Batch Method
        hi = Deconvolve_Batch(set_g);
    otherwise
        error('BOOL_DECONVOLVE must be either (Single) or (Batch)')
end




end
