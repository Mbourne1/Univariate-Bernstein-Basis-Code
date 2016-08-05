function hi = Deconvolve_Separate(set_g)
% Given a set of polynomials g_{i}, Perform a series of deconvolutions 
% independently, and output the solutions as an array of vectors.
%
% Inputs
%
%
% set_g : the array of polynomials on which deconvolution is performed.
%
%


% for each item in set g starting at
for i = 1:1:length(set_g)-1
    % Get the two polynomials which are to be deconvolved
    % f{i},f{i+1}
    hi{i} = Deconvolve(set_g{i},set_g{i+1});
end


end
