function [ t ] = o_roots_bisection( example_number,e_low,e_up, seed )
%   ROOTS_BISECTION obtain roots in interval by bisection method, once a
%   root is obtained with assumed multiplicity one, deconvolve and perform
%   bisection on f2 with root removed, until no more roots are found.

% bool_log
% 1 :   Use logs
% 0 :   Dont use logs.
global bool_log
bool_log = 1;

% Set Deconvolution Method
% 0 := Use standard non batch method
% 1 := Use batch deconvolution
global bool_deconvolve
bool_deconvolve = 0;

format long
eps_abs = 1e-10;
min_interval_size = 1e-10;
% Bool whether to print messages
bool_printMessageLog = 0;

% Get the exact roots
[f_roots_exact,~] = Root_Examples(example_number,seed);

% Get polynomial coefficients in scaled bernstein basis
f_exact_bi = B_poly(f_roots_exact);

% Get degree of fx
m = length(f_exact_bi) -1;

% Get binomial coefficients corresponding to f
Bi_n = zeros(m,1);
for i = 0:1:m
    Bi_n(i+1) = nchoosek(m,i);
end

% Obtain exact coefficients of polynomial f
f_exact = f_exact_bi./Bi_n;

% add noise to the coefficients
fx = VariableNoise(f_exact,e_low,e_up,seed);

% set interval
a = -1;
b = 1;


%while (b - a >= min_interval_size || ( abs( BernsteinEval(fx,a) ) >= eps_abs && abs( BernsteinEval(fx,b) )  >= eps_abs ) )
while (b - a >= min_interval_size && ( abs( BernsteinEval(fx,a) ) >= eps_abs && abs( BernsteinEval(fx,b) )  >= eps_abs ) )
    
    switch bool_printMessageLog
        case 1
            fprintf('Current Interval : %2.3f - %2.3f \n',a,b)
    end
    c = (a + b)/2;
    
    
    if ( abs(BernsteinEval(fx,c)) < eps_abs )
        switch bool_printMessageLog
            case 1
                fprintf('Root at %2.8f \n',c)
        end
        break;
    elseif ( BernsteinEval(fx,a)*BernsteinEval(fx,c) < 0 )
        switch bool_printMessageLog
            case 1
                fprintf('Change of sign in first half of bisection\n')
        end
        b = c;
    else
        a = c;
    end
end

if (b-a <= min_interval_size)
    fprintf('\nROOTS CALCULATED BY BISECTION FUNCTION \n');
    fprintf('No Roots Were Found\n')
    t = [];
    return
end

% Get calculated roots and multiplicity
r1 = [c 1];

% Initialise list of all roots
all_roots = r1;

% Convert to polynomial in scaled bernstein form
gx_bi = B_poly(r1);

% Get degree of polynomial formed from found roots
n = length(gx_bi) - 1;

% Calculate corresponding binomial coefficients
Bi_n = zeros(n+1,1);
for i=0:1:n
    Bi_n(i+1) = nchoosek(n,i);
end

% Build polynomial of removed root
gx = gx_bi ./ Bi_n;

% Perform deconvolution to obtain f2
deconvArray = {fx, gx};
f2 = Deconvolve(deconvArray);
f2 = cell2mat(f2(1));

while length(f2) ~=1
    a = 0;
    b = 1;
    
    
    while (b - a >= min_interval_size && ( abs( BernsteinEval(f2,a) ) >= eps_abs && abs( BernsteinEval(f2,b) )  >= eps_abs ) )
        switch bool_printMessageLog
            case 1
                fprintf('Current Interval : %2.3f - %2.3f \n',a,b)
        end
        c = (a + b)/2;
        if (  abs(BernsteinEval(f2,c)) < eps_abs )
            switch bool_printMessageLog
                case 1
                    fprintf('Found root : %2.8f\n', c)
            end
            break;
        elseif ( BernsteinEval(f2,a)*BernsteinEval(f2,c) < 0 )
            switch bool_printMessageLog
                case 1
                    
                    fprintf('Change of sign in first half of bisection\n')
            end
            b = c;
        elseif ( BernsteinEval(f2,b)*BernsteinEval(f2,c) < 0 )
            switch bool_printMessageLog
                case 1
                    fprintf('Change of sign in the second half of the bisection\n')
            end
            a = c;
        else
            switch bool_printMessageLog
                case 1
                    fprintf('Same sign in both, therefore - No roots found')
            end
            c = [];
            break
        end
    end
    
    % Get calculated roots and multiplicity
    
    r1 = [c 1];
    
    try
        % Add root to list of roots
        all_roots = [all_roots ;r1];
    catch
        break
    end
    
    % Convert to polynomial in scaled bernstein form
    gx_bi = B_poly(r1);
    
    % Get degree of polynomial formed from found roots
    n = length(gx_bi) - 1;
    
    % Calculate corresponding binomial coefficients
    Bi_n = zeros(n+1,1);
    for i=0:1:n
        Bi_n(i+1) = nchoosek(n,i);
    end
    
    % Build polynomial of removed root
    gx = gx_bi ./ Bi_n;
    
    % Perform deconvolution to obtain f2
    
    deconvArray = {f2, gx};
    f2 = Deconvolve(deconvArray);
    f2 = cell2mat(f2(1));
    
end

t = all_roots;

fprintf('\nROOTS CALCULATED BY BISECTION FUNCTION \n');
fprintf('\t Root \t \t \t \t\t \t \t \t \t \t \t \t Multiplicity \n')
fprintf('%22.15f \t \t\n',real(t(:,1)));
end

function [xx] =  BernsteinEval(f,c)
%% Evaluate function f at point c \in t
y = c;
m = length(f)-1;

for i = 0:length(f)-1
    
    x(i+1) = f(i+1) .* nchoosek(m,i) .* (1-y)^(m-i) .* y^i ;
end

xx = sum(x);

end


