function arr_hx = GetPolynomialArrayWithoutThetas(arr_hw, theta)
% Get array of polynomials h_{i}(x) from h_{i}(\omega)
%
% % Inputs
%
% arr_hw : (Array of Vectors)
%
% theta : (Float)
%
% % Outputs
%
% arr_hx : (Array of Vectors)

% Get number of polynomials in the array
nPolynomials_arr_hx = length(arr_hw);


% Initialise a cell array to store h_{i}(x)
arr_hx = cell(nPolynomials_arr_hx, 1);


for i = 1:1:nPolynomials_arr_hx
    
    arr_hx{i} = GetWithoutThetas(arr_hw{i}, theta);
    
end

end