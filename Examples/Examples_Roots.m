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



