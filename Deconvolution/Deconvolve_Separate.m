function arr_hx = Deconvolve_Separate(arr_fx)
% Given a set of polynomials g_{i}, Perform a series of deconvolutions 
% independently, and output the solutions as an array of vectors.
%
% Inputs
%
%
% arr_fx : the array of polynomials on which deconvolution is performed.
%
%

% Get number of polynomials in array of f_{i}(x)
nPolys_arr_fx = size(arr_fx,1);

% Initialise the array of h_{i}(x)
arr_hx = cell(nPolys_arr_fx - 1,1);

% for each item in set g starting at
for i = 1:1:nPolys_arr_fx - 1
    
    % Get the two polynomials which are to be deconvolved
    % f{i},f{i+1}
    arr_hx{i,1} = Deconvolve(arr_fx{i},arr_fx{i+1});
    
end


end
