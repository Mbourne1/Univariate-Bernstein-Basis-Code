
function [a,b,c,rankloss,u,v]=PreviousExamples(n,seed)

% This function is a database of pairs of polynomials. The pairs of
% polynomials are indexed by n.

% The matrices a and b define polynomials, where the first column of a
% and b defines the root, and the second column defines the multiplicity
% of the root.

% The matrix c stores the GCD of a and b, in the same format as a and b.

% Examples 1-20 Ning Examples
% Examples 20-40 My Examples
% Examples 100 - Report Examples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch n
    
    case 1
        a=[0.1,   3;
            0.5,   3;
            0.7,   3;
            1.3,   3;
            -1.7,   3;
            2.2,   3];
        
        b=[0.1,   4;
            0.5,   5;
            0.7,   2;
            1.0,   3;
            -1.0,   3;
            2.2,   5];
        
        
        
        
    case 2
        a=[
            -1.6,   5;
            0.1,   4;
            0.6,   5;
            0.9,   4;
            1.4,   5;
            2.0,   3];
        
        b=[
            -2.0,   3
            -1.6,   7;
            -1.4,   3;
            -0.1,   5;
            0.6,    4;
            0.9,    4;
            ];

        
        
    case 3
        a=[0.3,   5;
            -0.6,   4;
            -1.1,   5;
            2.7,   5;
            1.1,   2];
        
        b=[0.8,   6;
            -0.6,   5;
            -1.1,   3;
            2.7,   5;
            -0.8,   2];
       
        
        
    case 4
        a=[0.1,   4;
            0.5,   3;
            0.9,   4;
            1.4,   5;
            -0.7,   4;
            2.4,   3;
            1.0,   4];
        
        b=[0.1,   2;
            0.7,   2;
            0.9,   5;
            1.4,   5;
            -0.7,   4;
            -2.4,   3;
            -1.0,   4];

        
        
    case 5
        a=[0.14,  3;
            0.56,  3;
            0.89,  4;
            1.45,  4;
            2.37,  5;
            -3.61,  5];
        
        b=[0.14,  4;
            -0.76,  2;
            0.99,  1;
            -1.24,  2;
            2.37,  3;
            -3.61,  7];
        
        
        
        
    case 6
        a=[0.23,  5;
            0.56,  4;
            0.79,  3;
            -1.57,  5;
            3.61,  3;
            8.71,  5];
        
        b=[0.23,  6;
            0.56,  3;
            0.79,  2;
            -1.57,  6;
            -2.61,  2;
            8.71,  6];

        
        
        
    case 7
        a=[0.23,  5;
            0.56,  4;
            0.79,  3;
            -1.57,  3;
            -2.61,  3;
            4.71,  3;
            1.45,  4];
        
        b=[0.23,  4;
            -0.36,  3;
            -0.49,  2;
            1.32,  2;
            -2.61,  2;
            4.71,  3;
            1.45,  5];
        

        
        
    case 8
        a=[0.10,  3;
            0.56,  5;
            1.40,  4;
            -2.68,  3;
            1.79,  3;
            2.69,  2];
        
        b=[0.10,  4;
            0.56,  4;
            -1.40,  2;
            2.68,  2;
            1.79,  3;
            2.69,  3];
        
        
    case 9
        a=[0.10,  3;
            0.56,  4;
            0.75,  3;
            0.82,  3;
            1.37,  3;
            -0.27,  3;
            1.46,  2];
        
        b=[0.10,  2;
            0.56,  4;
            0.75,  3;
            0.99,  4;
            1.37,  3;
            2.12,  3;
            -1.20,  3];

        
        
    case 10
        a=[-0.20,  4;
            -0.60,  3;
            0.90,  4;
            9.60,  5;
            -4.30,  4;
            -2.00,  2];
        
        b=[0.20,   6;
            0.60,   2;
            0.90,   3;
            9.60,   5;
            -4.30,   3;
            2.60,   2];
        
        
    case 11
        a=[-0.10,  4;
            -0.50,  3;
            0.90,  4;
            9.60,  5;
            -4.30,  4;
            -2.00,  3];
        
        b=[0.20,   6;
            0.60,   2;
            0.90,   3;
            9.60,   5;
            -4.30,   3;
            2.60,   2];

        
    case 12
        a=[0.30,   6;
            0.60,   4;
            0.90,   6;
            9.70,   2;
            -4.60,   2;
            -2.30,   2];
        
        b=[0.30,   5;
            0.60,   5;
            0.90,   7;
            9.70,   2;
            4.60,   2;
            2.30,   2];
        

    case 13
        a=[0.23,   4;
            0.43,   3;
            0.57,   3;
            0.92,   3;
            1.70,   3];
        
        b=[0.23,   4;
            0.30,   2;
            0.77,   5;
            0.92,   2;
            1.20,   5];
        
        
    case 14
        a=[-0.10,  3;
            0.50,  4;
            0.70,  4;
            0.99,  3;
            -1.30,  3;
            1.50,  4];
        
        b=[-0.10,  4;
            0.50,  5;
            0.70,  2;
            0.90,  4;
            1.30,  2;
            1.50,  3];
        
        
    case 15
        a=[0.20,   3;
            0.60,   3;
            0.80,   3;
            -1.35,   2;
            1.90,   4];
        
        b=[0.20,   2;
            0.60,   5;
            0.80,   5;
            1.35,   3;
            -2.00,   3];
        
        
    case 16
        a=[-4.00,  2;
            0.10,  5;
            0.50,  3;
            1.40,  3;
            -1.60,  3];
        
        b=[0.70,   4;
            0.10,   3;
            -0.50,   2;
            1.40,   4;
            6.00,   4];
              
        
    case 17
        a=[0.10,   4;
            0.50,   5;
            0.90,   5;
            5.00,   6;
            -2.60,   4];
        
        b=[-0.10,  5;
            0.50,  4;
            0.90,  4;
            5.00,  4;
            2.00,  4];

    case 18
        a=[0.10,   4;
            0.50,   5;
            0.90,   5;
            5.00,   6;
            -2.60,   4];
        
        b=[-0.10,   5;
            -0.50,   1;
            0.90,   4;
            5.00,   4;
            2.00,   4];
               
        
    case 19
        a=[0.10,   3;
            0.40,   4;
            0.50,   4;
            0.70,   4;
            0.90,   3];
        
        b=[0.10,  5;
            0.40,  3;
            0.50,  6;
            0.70,  5;
            0.90,  4];
        
    % My examples     
        
    case 20
        a=[0.10,   3;
            0.20,   4;
            0.50,   4;
            0.70,   4;
            0.90,   3];
        
        b=[0.10,   5;
            0.40,   3;
            0.50,   6;
            0.70,   5;
            0.90,   4];

        
    case 21
        a = [
            0.10    7
            0.56    4
            0.40    4
            0.79    3
            0.69    2
            ];
        b = [
            0.10    9
            0.56    4
            0.69    2
            ];

        
        
    case 22
        a = [
            0.90    7
            0.56    4
            0.40    4
            0.79    3
            0.69    2
            ];
        b = [
            0.90    9
            0.56    4
            0.69    2
            ];

        
        
    case 23
        a = [
            1.10    7
            0.56    4
            0.40    4
            0.79    3
            0.69    2
            ];
        b = [
            1.10    9
            0.56    4
            0.69    2
            ];

        
        
    case 24
        a = [
            0.90    20
            0.60    10
            0.70    4
            0.40    35
            0.10    5
            ];
        b = [
            0.90    20
            0.60    10
            0.40    35
            0.10    10
            ];

        
        
    case 100 % Not a good example, but shows difference between build up and from scratch methods.
        a = [
            -2.68   3
            0.10    3
            0.56    4
            1.40    4
            1.79    3
            2.39    2
            ];
        
        b = [...
            -1.40   2
            0.10   4
            0.56   4
            1.79   3
            2.68   2
            
            ];
       
        
    case 101
        a = [
            0.14   3
            0.56   3
            0.89   4
            1.45   4
            2.37   3
            -3.61  5
            ];
        
        b = [
            0.14   4
            0.99   1
            2.37   3
            -0.76  2
            -1.24  2
            -3.61  7
            ];


        
    case 102
        a = [
            1.4630195735   6
            0.2120295267   5
            0.7231266078   5
            1.13713466634  4
            1.02953282453  4
            1.04908547995  5
            ];
        
        b = [
            1.4630195735    5
            0.2120295267    3
            5.134984201310244   4
            -1.659050925144660  5
            0.435719859785881   4
            0.759494024616640   3
            ];
    case 998
        a = [
            0.5 3
            0.7 2];
        b = [
            0.5 2
            0.7 1];
        
        
    case 999
        % relatively small example
        t = 3;
        m = 10;
        n = 7;
        [a,b] = BuildRandomPolynomials3(m,n,t,0,2,seed);
    
    case 1000
        t = 5;
        m = 15;
        n = 10;
        [a,b] = BuildRandomPolynomials3(m,n,t,0,2,seed);
        
    case 1001
        % increased degree of GCD from 1001
        t = 8;
        m = 15;
        n = 10;
        [a,b] = BuildRandomPolynomials3(m,n,t,0,2,seed);
        
    case 1002
        % increased degree of input polynomials m and n
        t = 10;
        m = 20;
        n = 15;
        [a,b] = BuildRandomPolynomials3(m,n,t,0,1,seed);
    
    
    % testing where deg(GCD) = deg(v) so that all sylvester subresultants
    % are rank deficient.
    
    case 1003;
        t = 4;
        m = 5;
        n = 4;
        [a,b] = BuildRandomPolynomials3(m,n,t,0,1,seed);
    
    case 1004
        t   =   7;
        m   =   10;
        n   =   9;
        [a,b] = BuildRandomPolynomials3(m,n,t,0,1,seed);
end

format long 




c = getDivisor(a,b);
u = getQuotient(a,c);
v = getQuotient(b,c);

% Sort roots of a into ascending order
[values, order] = sort(a(:,1));
a = a(order,:);

% Sort roots of b into ascending order
[values, order] = sort(b(:,1));
b = b(order,:);


m = sum(a(:,2));
n = sum(b(:,2));
d = sum(c(:,2));

rankloss = d;
rankloss = sum(c(:,2));

end