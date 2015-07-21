function [roots_calc] = o_roots_mymethod(ex_num,emin,emax,seed)
% o_roots_mymethod - Obtain polynomial roots of a noisy polynomial
% Given a set of polynomial roots, for a polynomial f, expressed in the
% Bernstein basis. A Bernstein basis polynomial f is produced, noise added,
% and its roots are returned, by method of finding GCDs of polynomials and
% their derivatives.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs

% ex - (Int) Example Number

% emin - Signal to noise ratio (minimum)

% emax - Signal to noise ratio (maximum)

% seed - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Global Variables

global output_format;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Add necessary paths
addpath 'Measures'
addpath 'Examples'
addpath 'multroot/multroot'

% Get coefficients in Scaled Bernstein Basis
[f_roots_exact,~]   = Root_Examples(ex_num,seed);
f_exact_bi          = B_poly(f_roots_exact);

% Get degree of f
m = length(f_exact_bi) - 1;

% Display the degree of the input polynomial
disp('Degree of Input Polynomial F ');
disp(int2str(m));

% Get the Binomial coefficients corresponding to the coefficients of
% polynomial f.
Bi_m = zeros(m+1,1);
for i=1:1:m+1
    Bi_m(i) = nchoosek(m,i-1);
end

% Get coefficients in Bernstein Basis
f_exact = f_exact_bi./Bi_m;

% Add Noise to coefficients of exact polynomial f_exact, to obtain noisy
% polynomial fx.
fx = VariableNoise(f_exact,emin,emax,seed);

% Obtain gx, the derivative of the noisy polynomial fx
gx = Differentiate_BernsteinBasis(fx)';


gx_exact = Differentiate_BernsteinBasis(f_exact)';

% Get gcd and quotients u and v, of fx and gx

% Initialise an iteration counter
ite_num = 1;

% Print current iteration number
fprintf('GCD Calculation Loop iteration = %i \n\n', ite_num);

% Perform an intial GCD Calculation, from this we obtain new forms of f and
% g and a calculated 'd'.
[f,g,d, ~ ,~,alpha,theta] = o1(fx,gx);


if output_format == 1
    % we have output in the form fw, gw, dw
    
    % divide GCD by theta vector to obtain dx
    dx = zeros(length(d)-1,1);
    for i = 1:1:length(d)
        dx(i) = d(i)./(theta^(i-1));
    end
    
    % divide f by theta vector to obtain fx
    fx = zeros(m+1,1);
    for i = 1:1:length(f)
        fx(i) = f(i)./(theta^(i-1));
    end
    
elseif output_format == 0
    % we have output in the form fx,gx,dx    
    fx = f;
    gx = g;
    dx = d;
end



% Initialise an array 'q' which stores the gcd outputs from each gcd
% calculation
q{1} = fx;
q{2} = dx;

% let degree_vector store the degrees corresponding to the array of
% GCDs stored in q.
degree_vec(1) = length(fx)-1;
degree_vec(2) = length(dx)-1;

% Let theta_vec store all theta values used in each iteration.
theta_vec(1) = 1;
theta_vec(2) = theta;

% increment the iteration number
ite_num = ite_num+1;


% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while length(q{ite_num})-1 > 0
    
    % set polynomial f to be the most recently calculated GCD.
    fx = q{ite_num};
    
    % set polynomial g to be the derivative of the GCD.
    gx = Differentiate_BernsteinBasis(q{ite_num})';
    
    % get degrees m and n of polynomials f and g respectively.
    m = size(fx,1) -1;
    n = size(gx,1) -1;
    
    % if degree of f is greater than one
    if m > 1
        fprintf('GCD Calculation Loop iteration = %i \n\n',ite_num );
        [f,g,d, u ,v,alpha,theta] = o1(fx,gx);

        
        if output_format == 1
            % we have output in the format fw, gw and dw
            
            % divide GCD by theta vector
            dx = zeros(length(d),1);
            for i = 1:1:length(d)
                dx(i) = d(i)./(theta^(i-1));
            end
            
            % divide f by theta vector
            fx = zeros(length(f),1);
            for i = 1:1:length(f)
                fx(i) = f(i)./(theta^(i-1));
            end
            
        elseif output_format == 0
            % we have output in the format fx, gx and dx
            
            fx = f;
            gx = g;
            dx = d;
            
        end
        
        % add the value of theta used in this GCD calculation to the theta
        % vector
        theta_vec(ite_num+1) = theta;
        
        % add the degree of the calculated GCD to the degree vector
        degree_vec(ite_num+1) = length(dx)-1;
        
        % replace input f with updated f
        q{ite_num} = fx;
        
        % add vector dx to the array of gcds 'q'
        q{ite_num+1} = dx;
        
        % increment iteration number.
        ite_num = ite_num+1;
        
        
    elseif m == 1
        % if m=1, then n = 0, GCD has maximum degree 0.
        fprintf('Only one subresultant exists \n')
        dx = 1;
        theta_vec(ite_num+1) = 1;
        degree_vec(ite_num+1) = length(dx)-1;
        q{ite_num+1} = dx;
        ite_num = ite_num+1;
        break;
        
        
    end
end

% #########################################################################
% This section performs deconvolution by structured matrix method.
addpath 'BernsteinMethods'

% fprintf('All Polynomials from GCD Calculations')
% for i = 1:1:length(q)
%     q{i}
% end

% Deconvolve the first set of polynomials.
h1 = Deconvolve(q);

% Deconvolve the second set of polynomials
w1 = Deconvolve(h1);


% w1 yields the simple, double, triple roots of input polynomial f.
% w1{i} yields the roots of multiplicity i.

% set the w1{max} = h1{max}
w1{ite_num-1} = h1{ite_num-1};

% get number of entries in w1
[~,c] = size(w1);

% initialise an empty set
wp1 = [];

% for each multiplicity in w1.
for i = 1:1:c
    
    % if the polynomial of said multiplicity is of length 2, degree 1, then
    % only one root exists for this multiplicity. Add it to the list wp1.
    if (length(w1{i}) == 2)
        wp = [];
        % get the polynomial, whose roots have multiplicty i, in bernstein form
        xx = w1{i}';
        % Normalise the polynomial
        xx = xx./xx(1);
        % Convert to power form
        yy = [xx(:,1) xx(:,2)-xx(:,1)];
        yy = [-yy(:,1)./yy(:,2) yy(:,2)./yy(:,2)];
        wp1 = [wp1 ; yy];
    elseif (length(w1{i}) > 2)
        % The given multiplicity contains more than one root, such that
        % number of coefficients in greater than 2, use MATLAB roots
        % function to find roots.
        % display('Multiplicity contains more than one root. ')
        
        w1{i}';
        n = length(w1{i}) -1;
        for y = 0:1:n
            bi_n(y+1) = nchoosek(n,y);
        end
        w1{i} = w1{i}'.* bi_n;
        xx = roots(w1{i});
        xxy = [1- xx./(1+xx) ];
        one = ones(length(xxy),1);
        xxy = [xxy one];
        wp1 = [wp1 ; xxy];
    end
end



% Obtaining multiplicities of the calculated roots
roots_multiplicty = [];
while sum(degree_vec) ~= 0
    % Get index of first zero.;
    index = min(find(degree_vec==0)) - 1;
    minus_vector = zeros(1,length(degree_vec));
    minus_vector(1,1:index) = index:-1:1;
    degree_vec = degree_vec - minus_vector;
    roots_multiplicty = [roots_multiplicty ; index];
end



roots_calc = [wp1(:,1) flipud(roots_multiplicty)];

%% Print the calculated roots and the corresponding multiplicities.
fprintf('\nROOTS CALCULATED BY MY METHOD \n');
fprintf('\t Root \t \t \t \t\t \t \t \t \t \t \t \t Multiplicity \n')
fprintf('%22.15f + %22.15f i  \t \t %3g \n',[real(roots_calc(:,1)),imag(roots_calc(:,1)),...
    roots_calc(:,2)]');
fprintf('\n');






end


