function [ t ] = o_roots_bisection(fx)
%   ROOTS_BISECTION obtain roots in interval by bisection method, once a
%   root is obtained with assumed multiplicity one, deconvolve and perform
%   bisection on f2 with root removed, until no more roots are found.



while (b - a >= min_interval_size && ( abs( Bernstein_Evaluate(fx,a) ) >= eps_abs && abs( Bernstein_Evaluate(fx,b) )  >= eps_abs ) )
    
    switch bool_printMessageLog
        case 1
            fprintf('Current Interval : %2.3f - %2.3f \n',a,b)
    end
    
    % Obtain the midpoint of the interval
    c = (a + b)/2;
    
    % check to see if the evaluation of the function f at the midpoint is 
    % within the margins of error of zero
    if ( abs(Bernstein_Evaluate(fx,c)) < eps_abs )
        switch bool_printMessageLog
            case 1
                fprintf('Root at %2.8f \n',c)
        end
        break;
    
    % Check for change of sign in the first half of the interval    
    elseif ( Bernstein_Evaluate(fx,a)*Bernstein_Evaluate(fx,c) < 0 )
              
        switch bool_printMessageLog
            case 1
                fprintf('Change of sign in first half of interval\n')
        end
        % if a change of sign exists, then root lies in the interval [a,c] and
        % must bisect the new interval. so set upper end of interval b = c.
        b = c;
        
    % Check for change of sign in the second half of the interval
    elseif ( Bernstein_Evaluate(fx,c)*Bernstein_Evaluate(fx,b) < 0)
        switch bool_printMessageLog
            case 1
                fprintf('Change of sign in second half of interval\n')
        end
        % if change of sign exists, then root lies in the interval [c,b]
        % and we must bisect the new interval. so set lower end of the
        % interval a equal to c.
        
        a = c;
       
    else
        switch bool_printMessageLog
            case 1
            fprintf('Same sign in both halves of the interval, therefore - No roots found')
        end
        fprintf('\nROOTS CALCULATED BY BISECTION FUNCTION \n');
        fprintf('No Roots Were Found\n')
        t = [];
    return
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

% Build polynomial of the removed root
gx = gx_bi ./ Bi_n;

% Perform deconvolution to obtain f2, the remainder of the polynomial now
% that the found root has been removed.
deconvArray = {fx, gx};
f2 = Deconvolve(deconvArray);
f2 = cell2mat(f2(1));

% while the degree of f2 is greater than 1
while length(f2) ~=1
    % Set interval lower bound
    a = intvl_lwr_bound;
    % Set interval upper bound
    b = intvl_uppr_bound;
    
    
    while (b - a >= min_interval_size && ( abs( Bernstein_Evaluate(f2,a) ) >= eps_abs && abs( Bernstein_Evaluate(f2,b) )  >= eps_abs ) )
        switch bool_printMessageLog
            case 1
                fprintf('Current Interval : %2.3f - %2.3f \n',a,b)
        end
        
        % Get the midpoint of the interval
        c = (a + b)/2;
        
        % Evaluate the function at the midpoint
        if (  abs(Bernstein_Evaluate(f2,c)) < eps_abs )
            switch bool_printMessageLog
                case 1
                    fprintf('Found root : %2.8f\n', c)
            end
            break;
            
        % Check for change of sign in first half of the interval
        elseif ( Bernstein_Evaluate(f2,a)*Bernstein_Evaluate(f2,c) < 0 )
            switch bool_printMessageLog
                case 1
                    
                    fprintf('Change of sign in first half of bisection\n')
            end
            % if change of sign, then set the upper bound of the next
            % interval b = c
            b = c;
            
        % Check for change of sign in the second half of the interval
        elseif ( Bernstein_Evaluate(f2,b)*Bernstein_Evaluate(f2,c) < 0 )
            switch bool_printMessageLog
                case 1
                    fprintf('Change of sign in the second half of the bisection\n')
            end
            % if change of sign, then set the lower bound of the next
            % interval a = c
            a = c;
        else
            switch bool_printMessageLog
                case 1
                    fprintf('Same sign in both halves of the interval, therefore - No roots found')
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

% Print out roots
PrintoutRoots('BISECTION' , t)
end









