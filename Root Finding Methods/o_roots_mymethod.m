function [roots_calc] = o_roots_mymethod(fx)
% Calculate the roots of the input polynomial f(x) in Bernstein form.
%
% Inputs.
%
% fx:   Vector of polynomial coefficients in Bernstein basis, where the
%       ith entry of fx, fx(i) is the coefficient a_{i}B^{m}_{i}, and m is the
%       degree of fx.
%
% Outputs.
%
% roots_calc : The Calculated roots of the polynomial f(x)


addpath 'Measures'
addpath 'Root Finding Methods/multroot/multroot'


% Initialise an iteration counter
ite_num = 1;

% Initialise an array 'q' which stores the gcd outputs from each gcd
% calculation
q{1} = fx;

% let degree_vector store the degrees corresponding to the array of
% GCDs stored in q.
degree_vec(1) = size(fx,1)-1;

% Let theta_vec store all theta values used in each iteration.
theta_vec(1) = 1;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while length(q{ite_num})-1 > 0
    
    % set polynomial f to be the most recently calculated GCD.
    fx = q{ite_num};
    
    % set polynomial g to be the derivative of the GCD.
    gx = Bernstein_Differentiate(q{ite_num});
    
    % get degrees m and n of polynomials f and g respectively.
    m = GetDegree(fx);
    
    fprintf('Begin : GCD Calculation Loop iteration = %i \n\n',ite_num );
    fprintf('m = %i \n',m);
    fprintf('n = %i \n\n',m-1);
    
    % if degree of f is greater than one
    if m > 1
        
        
        [fx,~,dx, ~ ,~,~,theta] = o1(fx,gx);
        
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
        %theta_vec(ite_num+1) = 1;
        degree_vec(ite_num+1) = length(dx)-1;
        q{ite_num+1} = dx;
        ite_num = ite_num+1;
        break;
        
        
    end
end

% #########################################################################
% This section performs deconvolution by structured matrix method.

% fprintf('All Polynomials from GCD Calculations')
% for i = 1:1:length(q)
%     q{i}
% end

% Deconvolve the first set of polynomials.
h1 = Deconvolve(q);

[~,nEntries_h] = size(h1);
if nEntries_h == 1
    
    % if number of cols in h1 is only 1, then do not perform second set of
    % deconvolutions, since only one entry in h1.
    
    % Initialise an empty root array
    root_arr = [];
    
    % Get the polynomial a(w)
    aw = h1{1};
    
    % Normalise the coefficients of a(w)
    aw = aw./aw(1);
    
    % Get the degree of a(w)
    deg_aw = GetDegree(h1{1});
    
    % get aw including its binomial coefficients.
    aw_bi = GetWithBinomials(aw);
    
    % get the roots in terms of z^{i} where z^{i} =(\frac{y}{(1-y)})^{i}.
    % Note we must flip the coefficients to conform with input format of
    % MATLAB roots function.
    rt_wrt_z = roots(flipud(aw_bi));
    
    % Get the roots in terms of y.
    roots_wrt_y = [rt_wrt_z./(1.+rt_wrt_z) ]; %Edit 27/07
    
    % Initialise a vector of ones
    one = ones(length(roots_wrt_y),1);
    
    % Get the roots with respect to y, and their multiplicities all set
    % to one.
    roots_wrt_y = [roots_wrt_y one];
        
    % Add the roots to the array of roots
    root_arr = [root_arr ; roots_wrt_y];
else
    
    % perform deconvolutions
    
    % Deconvolve the second set of polynomials h_{i} to obtain the set of
    % polynomials w_{i}
    w1 = Deconvolve(h1);
    
    
    % w1 yields the simple, double, triple roots of input polynomial f.
    % w1{i} yields the roots of multiplicity i.
    
    % set the w1{max} = h1{max}
    w1{ite_num-1} = h1{ite_num-1};
    
    % get number of entries in w1
    [~,nEntries_w1] = size(w1);
    
    % initialise an empty set
    root_arr = [];
    
    % for each multiplicity in w1.
    for i = 1:1:nEntries_w1
        
        
        % if the polynomial of said multiplicity is of length 2, degree 1, then
        % only one root exists for this multiplicity. Add it to the list wp1.
        if (length(w1{i}) == 2)
            
            
            % Get the polynomial a(w), whose simple roots have multiplicty 
            % i in the polynomial f which we started with.
            aw = w1{i};
            
            % Normalise the polynomial coefficients
            aw = aw./aw(1);
            
            % Convert to power form, so that coefficients are in terms of y^{i}
            % rather than (1-y)^{m-i}y^{i}.
            a_pwr = [aw(1,:) ; aw(2,:)-aw(1,:)];
            
            % Obtain the root in terms of y, and set multiplicity to one.
            a_rt = [-a_pwr(1,:)./a_pwr(2,:) a_pwr(2,:)./a_pwr(2,:)];
            
            % Add the root to the [root, mult] matrix
            root_arr = [root_arr ; a_rt];
            
        elseif (length(w1{i}) > 2)
            % The given multiplicity contains more than one root, such that
            % number of coefficients in greater than 2, use MATLAB roots
            % function to find roots.
            % display('Multiplicity contains more than one root. ')
            
            % get the polynomial a(w), whose roots have multiplicity i, in bernstein
            % form.
            aw = w1{i};
            
            % Normalise the polynomial coefficients
            aw = aw./aw(1);
            
            % get the degree of a(w)
            n = GetDegree(w1{i});
            
            % get aw including its binomial coefficients.
            aw_bi = GetWithBinomials(aw);
            
            % get the roots in terms of z^{i} where z^{i} =(\frac{y}{(1-y)})^{i}.
            rt_wrt_z = roots(flipud(aw_bi));
            
            % get the root in terms of y
            roots_wrt_y = [rt_wrt_z./(1.+rt_wrt_z) ]; %Edit 27/07
            
            % Initialise a vector of ones
            one = ones(length(roots_wrt_y),1);
            
            % get the roots with respect to y, and their multiplicities all set
            % to one.
            roots_wrt_y = [roots_wrt_y one];
            
            % add the roots to the array of roots
            root_arr = [root_arr ; roots_wrt_y];
        end
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



roots_calc = [root_arr(:,1) flipud(roots_multiplicty)];

%% Print the calculated roots and the corresponding multiplicities.
PrintoutRoots('MY METHOD',roots_calc);



