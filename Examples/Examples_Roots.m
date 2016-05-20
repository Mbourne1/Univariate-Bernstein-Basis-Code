function [fx_exact] = Examples_Roots(ex_num)

global Example_Type

switch Example_Type
    case 'From Roots'
        
        
        % Get exact polynomial roots from the example file.
        [f_roots_exact,~] = Examples_Roots_FromRoots(ex_num);
        f_roots_exact = sortrows(f_roots_exact,1);
        
        % Print the exact roots and multiplicities to terminal
        fprintf('\nExact Roots of Input Polynomial \n');
        fprintf('\t \t \t \t \t Root \t \t \t \t\t  Multiplicity \n')
        fprintf('%30.15f %30.15f \t \t\n',[f_roots_exact(:,1),f_roots_exact(:,2)]');
        fprintf('\n');
        
        fx_bi = B_poly(f_roots_exact);
        fx_exact = GetWithoutBinomials(fx_bi);
        
        % Get degree of f
        m = GetDegree(fx_bi);
        
        % Display the degree of the input polynomial
        disp('Degree of Input Polynomial F ');
        disp(int2str(m));
        
        
    case 'From Coefficients'
        fx_exact = Examples_Roots_FromCoefficients(ex_num);
        %fx_bi = GetWithBinomials(fx_exact);
    otherwise
        error('error : Example_Type not valid')
end
end


function [root_mult_array_fx,m] = Examples_Roots_FromRoots(ex_num)
% Examples_Roots(n)
%
% Get a set of roots and multiplicities of a polynomial f(x).
%
% Inputs.
%
%
% n : Example Number
%
% Outputs
%
% root_mult_array_fx : Matrix of roots of f(x) and their corresponding
%                      multiplicities. Each row is a root, multiplicity
%                      pair.
%
% m : Degree of polynomial f(x).


switch ex_num
    
    case 'Example'
        root_mult_array_fx =...
            [
            0.1     7;
            0.9     12;
            ];
        
    case 'Example Zeng'
        root_mult_array_fx = ...
            [
            10/11   5
            20/11   3
            30/11   2
            ];
        
    case '-1'
        root_mult_array_fx = ...
            [
            0.1     1;
            0.5     10;
            ];
        
    case '0'
        root_mult_array_fx = ...
            [
            0.1     1;
            0.7     1;
            0.9     1;
            ];
        
    case '1'
        root_mult_array_fx = ...
            [
            0.1     1;
            0.9     2;
            0.2     3;
            0.5     4;
            0.6     5;
            ];
        
    case '2'
        root_mult_array_fx = ...
            [
            0.5     1;
            0.8     2;
            ];
        
    case '3'
        root_mult_array_fx = ...
            [
            0.3,     1;
            0.6,    2;
            -1.1,    3;
            ];
    case '4'
        root_mult_array_fx = ...
            [
            0.1    1;
            0.5    2;
            0.9    3;
            1.4    4;
            ];
        
    case '5'
        root_mult_array_fx = ...
            [
            0.1   1;
            0.5   2;
            0.8   3;
            0.4   4;
            1.3   5;
            ];
        
    case '6'
        root_mult_array_fx = ...
            [
            0.14,  1;
            0.56,  2;
            0.89,  3;
            0.45,  4;
            0.37,  5;
            1.3,  6
            ];
        
    case '7'
        root_mult_array_fx = ...
            [
            0.14    1;
            0.56    2;
            0.89    3;
            1.45    4;
            2.37    5;
            -3.6    6;
            0.8     7;
            ];
        
    case '8'
        root_mult_array_fx = ...
            [
            0.14,  1;
            0.56,  2;
            0.89,  3;
            1.45,  4;
            2.37,  5;
            -3.61,  6
            0.8     7
            0.6     8
            ];
        
    case '9'
        root_mult_array_fx= ...
            [
            0.14,  1;
            0.56,  2;
            0.89,  3;
            1.45,  4;
            2.37,  5;
            -3.61,  6
            0.8     7
            0.6     8
            1.2     9
            ];
    case '10'
        root_mult_array_fx = ...
            [
            0.14,  1;
            0.56,  2;
            0.89,  3;
            1.45,  4;
            2.37,  5;
            -3.61,  6
            0.8     7
            0.6     8
            1.2     9
            0.10003245657   10
            ];
        
        % cases with roots of same multiplicites
    case '11'
        root_mult_array_fx = ...
            [
            0.1     1;
            0.4     2;
            0.7     2;
            ];
        
    case '12'
        root_mult_array_fx = ...
            [
            0.1     1;
            0.5     2;
            0.6     2;
            0.9     3;
            1.1     3;
            ];
        
    case '13'
        root_mult_array_fx = ...
            [
            0.1     1;
            0.5     2;
            0.6     2;
            0.9     3;
            1.1     5;
            ];
        
    case '14'
        root_mult_array_fx = ...
            [
            0.1     1;
            0.5     2;
            0.6     2;
            0.9     3;
            1.1     7;
            ];
        
    case 'Custom'
        
        prompt = 'Enter the degree of Polynomial f(x) :';
        m = input(prompt);
        
        intvl_low = -1;
        intvl_high = +1;
        
        root_mult_array_fx = BuildRandomPolynomial(m,intvl_low,intvl_high);
    otherwise
        error('Example Number not found')
        
end


m = sum(root_mult_array_fx(:,2)-1);

end

function fx_exact = Examples_Roots_FromCoefficients(ex_num)
switch ex_num
    case '1'
        fx_exact = ...
            [
            -0.9865
            2.2398
            2.8950
            1.9092
            -0.1477
            ];
    otherwise
        error('Not a valid example number for the *from coefficients* examples.')
end
end