function [gm]=GeometricMean(fx,n,k)
% This function calculates the geometric mean gm of the terms that
% contain the coefficients of the polynomial c in the kth subresultant
% matrix S(c,d). The integer n is the degree of the polynomial d.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%                            Inputs

% fx    :  The coefficients of one of the Bernstein basis polynomials,
%            f or g.

% n     :  The degree of the other polynomial.

% k     :  The order of the subresultant matrix.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                          Outputs


% gm :  The geometric mean of the terms in the modified Sylvester
%            Sylvester matrix S(f,g)Q.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Variables.

% BOOL_LOG - (Boolean)
%   1 :- Perform calculations by log method
%   0 :- Perform calculations by standard method.

global bool_log

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%S

% Get previous state of logs
previous_log_state = bool_log;

% To avoid overflow always set BOOL_LOG to be 1.
bool_log = 1;


% Dependent on which method is used, 1 - use logs, 0 - use nchoosek

switch bool_log
    case 1 % Calculate GM using logs
        gm = GMlog(fx,n,k);       
    case 0 % Calculate GM without logs
        gm = GMnchoosek(fx,n,k);
end

% reset log method
bool_log = previous_log_state;

end


function gm = GMlog(fx,n,k)
% Calculate the Geometric mean of the entries of the Coefficient matrix $C_{f}$
% which may or may not contain Q. may or may not contain denominator.#

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs.


% fx :  Polynomial coefficients of fx

% n :   Degree of polynomial gx

% k :   Index of subresultants S_{k}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Outputs.


% gm :  Geometric mean of entries of fx in the Syvlester Matrix S_{k}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Variables.


global bool_denom_syl

global bool_q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch bool_q
    case 1 
        % Include Q in Sylvester Subresultant S_{k} so include Q in the 
        % calculation of the geometric mean.
        
        % Calculate the degree of the polynomial.
        m = length(fx)-1;
        
        % Calculate the absolute value of the coefficients in c.
        fx=abs(fx);
        
        % Calculate part 1 of the geometric mean, the coefficient part.
        p1_log = (1/(m+1)).*log10(sum(fx));
        %p1_exp = 10^p1_log;
        
        % Calculate part 2 of the geometric mean.
        p2_log = 0;
        for j = 0:1:n-k
            for i = 0:1:m
                p2_log = p2_log + lnnchoosek(i+j,j);
            end
        end
        
        p2_log = (2./((n-k+1)*(m+1)))* (p2_log);
        %p2_exp = 10^p2_log;
        
        switch bool_denom_syl
            case 1 % if denominator is included
                %p3 = nchoosek(m+n-k,n-k);
                p3_log = lnnchoosek(m+n-k,n-k);
                %p3_exp = 10.^p3_log;
            case 0 % if denominator is not included
                p3_log = 0;
        end

        gm = 10.^(p1_log + p2_log - p3_log);

    case 0
        % Exclude Q from Geometric mean calculations.
        % Split this calculation in to three parts, Numerator_A, the
        % coefficients of fx, Numerator_B : the binomial coefficients
        % corresponding to the a_{i} in the scaled bernstein basis, and
        % Denominator. When Sylvester Matrix is without Q, numerators
        % consist of a_{i}\binom{m}{i}. Denominator changes for each row
        % and is given by \binom{m+n-k}{i+j} where i is the row number and
        % j is the column number, (index from 0).
        
        % Geometric mean is given by 
        % $$\prod_{j=0}^{n-k}\prod_{i=0}^{m}$$
        %
        %
        
        % Calculate the degree of the polynomial.
        m = length(fx)-1;
        
        % Calculate the absolute value of the coefficients in c.
        fx=abs(fx);
        
        % First Dealing with the numerators
        
        Numerator_A = 1;
        Numerator_B = 0;
        
        % For each coefficient a_{i} 
        for i = 0:1:m
            Numerator_A = Numerator_A .* fx(i+1);
            Numerator_B = Numerator_B + lnnchoosek(m,i);
        end
        
        GM_Numerator_A = Numerator_A .^(1/(m+1));
        GM_Numerator_B = (1/(m+1)) * Numerator_B;
        GM_Numerator_B = 10.^ GM_Numerator_B;
        GM_Numerator = GM_Numerator_A * GM_Numerator_B;
        
        
        
        Denominator_LOG = 0;
        for i = 0:1:m
            for j = 0:1:n-k
                Denominator_LOG = Denominator_LOG + lnnchoosek(m+n-k,i+j);
            end
        end
        
        
        GM_Denominator_LOG = (1/((m+1)*(n-k+1))) * Denominator_LOG;
        
        GM_Denominator = 10.^GM_Denominator_LOG;
        
        gm = GM_Numerator./GM_Denominator;
        
        
        
end

end


function gm = GMnchoosek(fx,n,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs.


% fx :  Polynomial coefficients of fx

% n :   Degree of polynomial gx

% k :   Index of subresultants S_{k}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Outputs.


% gm :  Geometric mean of entries of fx in the Syvlester Matrix S_{k}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Variables.


global bool_denom_syl

global bool_q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch bool_q
    
    case 1 % Include Q
        
        % Calculate the degree of the polynomial.
        m = length(fx)-1;
        
        p2 = 1;
        
        p1 = 1;
        for i = 0:1:m
            p1 = p1 .* (fx(i+1)^(1./(m+1)));
        end
        
        
        % since the product of the binomial coefficient A in the numerator
        % is equal to the product of the binomial coefficient B in the
        % numerator, only calulate this once.
        for j = 0:1:n-k
            for i = 0:1:m
                p2 = p2 .* (nchoosek(i+j,j)^(2./((n-k+1)*(m+1))));
            end
        end
 
        switch bool_denom_syl
            case 1 % if denominator is included
                p3 = nchoosek(m+n-k,n-k);
            case 0 % if denominator is not included
                p3 = 1;
        end
        
        gm = p1.*p2 ./ p3;
        
        
    case 0 % exclude Q
        % Calculate the degree of the polynomial f.
        m = length(fx)-1;
        
        % Initialise the product of the numerator at 1
        prod_num = 1;
        
        for i=0:1:m
            prod_num = prod_num .* fx(i+1) .* nchoosek(m,i);
        end
        
        % Initialise the product of the denominator at 1.
        prod_denom =1 ;
        
        for i = 0:1:m
            for j = 0:1:n-k        
                prod_denom = prod_denom * nchoosek(m+n-k,i+j);
            end
        end
        
        num =prod_num .^(1./(m+1));
        denom = prod_denom .^(1./((m+1)*(n-k+1)));
        
        gm = num/denom;
end
end

