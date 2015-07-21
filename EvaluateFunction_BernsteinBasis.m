function [f] = EvaluateFunction_BernsteinBasis(a,b,inc,fx)
% Evaluate the function f(x) over interval [a,b].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs.

% a :   lower limit of interval

% b :   upper limit of interval

% inc : size of steps within the interval [a,b]

% fx :  Coefficients of Polynomial f in the bernstein basis.

% f :   Vector of values corresponding to the evaluated function f, at 
% points x = [a:inc:b]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 0;

x = a:inc:b;
f = zeros(1,length(x));

for c = a:inc:b
    i = i+1;
    f(i) = Bernstein_Eval(fx,c);
end

end

function [xx] =  Bernstein_Eval(f,c)
% Evaluate function f at point c.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs.


% f :   Polynomial f

% c :   Point at which we evaluate polynomial f.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 y = c;
 m = length(f)-1;
 x = zeros(length(f),1);
 
    for i = 0:length(f)-1
        x(i+1) = f(i+1) .* nchoosek(m,i) .* (1-y)^(m-i) .* y^i ;
    end
    
    xx = sum(x);
    
end