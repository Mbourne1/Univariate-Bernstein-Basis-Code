

function [a,t]=RootsExamples(n,seed)
    
    % This function is a database of pairs of polynomials. The pairs of
    % polynomials are indexed by n.
    
    % The matrices a and b define polynomials, where the first column of a
    % and b defines the root, and the second column defines the multiplicity
    % of the root.
    
    % The matrix c stores the GCD of a and b, in the same format as a and b.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch n
        
        case -1
            a = [ 
                0.1     1;
                0.5     1;
                ];
            
        case 0
            a = [
                0.1     1;
                0.7     2;
                ];
            
        case 1
            a = [
                0.1     1;
                0.9     2;
                0.2     3;
                0.5     4;
                0.6     5;
                ];
            
        case 2
            a = [
                0.5     1;
                0.8     2;
                ];
            
        case 3
            a=[0.3,     1;
               0.6,    2;
               -1.1,    3;
               ];
        case 4
            a=[
                0.1    1;
                0.5    2;
                0.9    3;
                1.4    4;
                ];
            
        case 5
            a=[
                0.1   1;
                0.5   2;
                0.8   3;
                0.4   4;
                1.3   5;
                ];
            
        case 6
            a=[0.14,  1;
                0.56,  2;
                0.89,  3;
                0.45,  4;
                0.37,  5;
                1.3,  6
                ];
            
        case 7
            a=[
                0.14    1;
                0.56    2;
                0.89    3;
                1.45    4;
                2.37    5;
                -3.6    6;
                0.8     7;
                ];
            
        case 8
            a=[0.14,  1;
                0.56,  2;
                0.89,  3;
                1.45,  4;
                2.37,  5;
                -3.61,  6
                0.8     7
                0.6     8
                ];
            
        case 9
            a=[0.14,  1;
                0.56,  2;
                0.89,  3;
                1.45,  4;
                2.37,  5;
                -3.61,  6
                0.8     7
                0.6     8
                1.2     9
                ];
        case 10
            a=[0.14,  1;
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
        case 11
            a = [
                0.1     1;
                0.4     2;
                0.7     2;
                ];
            
        case 12
            a = [
                0.1     1;
                0.5     2;
                0.6     2;
                0.9     3;
            
            ];
            
        case 20
            m = 4;
            intvl_low = 0;
            intvl_high = 1;
            a = BuildPolynomial(m,intvl_low, intvl_high,seed);
            
         case 21
            m = 6;
            intvl_low = -1;
            intvl_high = 1;
            a = BuildPolynomial(m,intvl_low, intvl_high,seed);
            
         case 22
            m = 10;
            intvl_low = 0;
            intvl_high = 1;
            a = BuildPolynomial(m,intvl_low, intvl_high,seed);
            
         case 23
            m = 15;
            intvl_low = -1;
            intvl_high = 1;
            a = BuildPolynomial(m,intvl_low, intvl_high,seed);
            
         case 24
            m = 20;
            intvl_low = -1;
            intvl_high = 10;
            a = BuildPolynomial(m,intvl_low, intvl_high,seed);   
    end
    
    
    t = sum(a(:,2)-1);
    
end

