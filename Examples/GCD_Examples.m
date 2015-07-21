
function [a,b,c,rankloss,u,v]=GCD_Examples(n)
% Inputs.
% n - Index of example to be used
% Outputs.
% a - Roots and multiplicities of polynomial f.
% b - Roots and multiplicities of polynomial g.
% c - Roots and multiplicities of polynomial d, the GCD of f and g.
% u - Roots and multiplicities of quotient polynomial f/u = d.
% v - Roots and multiplicities of quotient polynomial g/v = d.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch n
    
    case -2 % From Bini
        a = [
            -1      1;
            -2     1;
            2      1;
            1      1;
            ];
        
        b = [
            2      1;
            -0.5   1;
            0.5    1;
            ];
    
    case -1
        a = [
            -1      1;
            -2     1;
            2      1;
            1      1;
            3      1;
            -3     1;
            ];
        
        b = [
            2      1;
            -0.5   1;
            0.5    1;
            3      1;
            -3     1;
            ];

    case 0
        a = [
            0.10    30
            0.56    4
            0.40    4
            0.79    3
            0.69    2
            ];
        
        b = [
            0.10    40
            0.56    4
            0.69    2
            ];

        
        
    case 1
        a = [
            0.2 2
            0.5 1
            0.7 1
            0.9 1
            1.1 1
            ];
        b = [
            0.5 1
            0.1 1
            0.9 1
            0.3 1
            ];

        
    case 2
        a = [
            0.2 1
            0.4 1
            0.6 1
            0.8 1
            ];
        
        b = [0.2 1
            0.3 1];
               
    case 3
        a = [
            0.56 20
            0.75 3
            0.82 3
            0.37 3];
        
        b = [
            0.56    20
            0.75    3
            0.99    4
            0.37    3
            0.12    3
            0.20    3
            ];
     
    case 4
        
        a = [0.1    20
            0.5    2];
        
        b = [0.1    20
            0.9    1];
        
    case 5
        
        a = [0.1    2
            0.3     2
            0.5    2];
        
        b = [0.1    2];
        
        % From Winkler Paper - Methods for the computation of the degree of an
        % approximate greatest common divisor of two inexact bernstein basis
        % polynomials.
    case 6
        a = [
            0.10    3
            0.56    4
            0.75    3
            0.82    3
            1.37    3
            -0.27   3
            1.46    2
            ];
        
        b = [
            0.10    2
            0.56    4
            0.75    3
            0.99    4
            1.37    3
            2.12    3
            1.20    3
            ];
              
    case 7
        a = [
            0.23   4
            0.43   3
            0.57   3
            0.92   3
            1.70   3
            ];
        b = [
            0.23   4
            0.30   2
            0.77   5
            0.92   2
            1.20   5
            ];
         
    case 8
        a = [0.1  5
            0.56 4
            0.75 3
            0.82 3
            1.37 3];
        b = [0.1  5
            0.56 4
            0.75 3
            0.99 4
            1.37 3
            2.12 3
            1.20 3];

    % Example 6.2
    case 9
        a = [
            0.14    3
            0.56    3
            0.89    4
            1.45    4
            2.37    3
            -3.61   5
            ];
        b = [
            0.14    4
            0.99    1
            2.37    3
            -0.76   2
            -1.24   2
            -3.61   7
            ];

    case 10
        a = [
            0.14    3
            0.56    3
            0.89    4
            0.37    3
            ];
        
        b = [
            0.14    3
            0.37    3
            0.76   2
            0.24   2
            ];
        
        
    % Example 6.2
    case 11
        a = [
            0.10    3
            0.56    5
            1.40    4
            1.79    3
            2.69    2
            -2.68   3
            ];
        b = [
            0.10    4
            0.56    4
            1.79    3
            2.68    2
            2.69    3
            -1.40   2
            ];
       
        
    case 999
        t = 5;
        m = 10;
        n = 7;
        [a,b] = BuildRandomPolynomials(m,n,t)
        
        
end
c = getDivisor(a,b)
u = getQuotient(a,c);
v = getQuotient(b,c);

m = sum(a(:,2));
n = sum(b(:,2));
d = sum(c(:,2));
rankloss = d;

end




