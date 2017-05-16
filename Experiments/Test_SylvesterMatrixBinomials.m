function [] = Test_SylvesterMatrixBinomials()
% This test looks at the sum of the diagonal entries of the Sylvester
% subresultant matrices.
%
%

% % Inputs
%
% m : (Int) Degree of polynomial f(x)
%
% n_k : (Int) Degree of polynomial v(x)

max_m = 5;
max_n_k = 5;

matrix_sumDiagonals(max_m + 1, max_n_k + 1);

for m = 0:1:max_m
    for n_k = 0:1:max_n_k
        
        Tf = zeros(m + n_k + 1, n_k+1);
        
        matrix_sumDiagonals(m+1,n_k+1) = test2(m,n_k);
        
        
        
    end
    
end

display(matrix_sumDiagonals);

end


function [my_sum] = test2(m,n_k)
%
% % Inputs
%
% m : (Int)
%
% n_k : (Int)


% Initialise matrix
my_vector = zeros(m+1,1);

for i = 0:1:m
    
    temp_sum = 0;
    
    for j = 0:1:n_k
        
%         coefficient_multiplier = ...
%             nchoosek(m,i) ...
%             * nchoosek(n_k,j) ...
%             * nchoosek(m+n_k,i+j);
        
         coefficient_multiplier = ...
             nchoosek(i+j,i)...
             * nchoosek(m+n_k-i-j,m-i)...
             / nchoosek(m+n_k,n_k);
        
        temp_sum = temp_sum + coefficient_multiplier;
    
        
        
    end
    
    my_vector(i+1) = temp_sum;
    
end

display(my_vector)
my_sum = temp_sum;

test_sum = (m+n_k+1) / (m+1) ;
display(my_sum)
display(test_sum)
end


