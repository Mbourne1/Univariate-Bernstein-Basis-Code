function [root_mult_array] = o_roots_mymethod(fx)
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
% root_mult_array : The Calculated roots of the polynomial f(x) and their
% corresponding multiplicities.

% Initialise an iteration counter
ite = 1;

% Initialise an array 'q' which stores the gcd outputs from each gcd
% calculation
f{1} = fx;

% let degree_vector store the degrees corresponding to the array of
% GCDs stored in q.
M(1) = GetDegree(f{1});

% Let theta_vec store all theta values used in each iteration.
vTheta(1) = 1;

% Get the number of distinct roots of f_{1}. Since this is unknown at this
% time, set number of distinct roots to be m_{1} = deg(f_{1}).
nDistinctRoots(1) = GetDegree(f{1});


% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while M(ite) > 0
    
   
    % if degree of f_{i} is greater than one
    if M(ite) > 1
        
        
        fprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite);
        
        % Get upper and lower bounds of the GCD Computation.
        % M_{i+1} > M_{i} - d_{i-1}
        try
            lower_lim = M(ite)-nDistinctRoots(ite-1);
            upper_lim = M(ite)-1;
            fprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim);
            fprintf('Maximum degree of f_{%i}: %i \n\n', ite+1, upper_lim);
        catch
            lower_lim = 1;
            upper_lim = M(ite)-1;
        end
    
        
        bool_canbe_coprime = false;
        
        % Get GCD of f(x) and f'(x)
        [f{ite},~,f{ite+1}, h1{ite} ,~,~,vTheta(ite+1), M(ite+1)] = o_gcd_mymethod( f{ite},...
            Bernstein_Differentiate(f{ite}), [lower_lim, upper_lim],bool_canbe_coprime);
        
        % Get number of distinct roots of f(ite)
        nDistinctRoots(ite) = M(ite) - M(ite+1);
        
        fprintf('Degree of f_{%i} : %i \n',ite + 1, M(ite+1))
        fprintf('Number of distinct roots in f_{%i} : %i \n',ite,nDistinctRoots(ite))

        % increment iteration number.
        ite = ite+1;
        
        
    elseif M(ite) == 1
        % if m=1, then n = 0, GCD has maximum degree 0.
        fprintf([mfilename ' : ' 'Only one subresultant exists \n'])
        dx = 1;
        %theta_vec(ite_num+1) = 1;
        M(ite+1) = length(dx)-1;
        f{ite+1} = dx;
        h1{ite} = f{ite};
        ite = ite+1;
        break;
        
        
    end
end


% Get the degree structure of the polynomials h_{i}
deg_struct_h = diff([M]);

% Get the degree structure of the polynomials w_{i}
deg_struct_w = diff([deg_struct_h 0]);

vMultiplicities = find(deg_struct_w~=0);


% #########################################################################
% This section performs deconvolution by structured matrix method.

% fprintf('All Polynomials from GCD Calculations')
% for i = 1:1:length(q)
%     q{i}
% end

% Deconvolve the first set of polynomials.

global SETTINGS
switch SETTINGS.ROOTS_HX
    case 'From Deconvolution'
        h1 = Deconvolve_Set(f);
    case 'From ux'
        h1 = h1;
end

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
    w1 = Deconvolve_Set(h1);
    
    
    % w1 yields the simple, double, triple roots of input polynomial f.
    % w1{i} yields the roots of multiplicity i.
    
    % set the w1{max} = h1{max}
    w1{ite-1} = h1{ite-1};
    
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
% create a matrix where the first column contains the multiplicities, and
% the second column contains the number of roots of that multiplicity
nPolys_wi = length(deg_struct_w);
mat = [(1:1:nPolys_wi)' deg_struct_w'];

count = 1;
root_mult_array = [];
for i = 1:1:size(mat,1)
    % Get the number of roots of multiplicity i
    nRoots_of_Multplicity_i = mat(i,2);
    
    for j = 1:1:nRoots_of_Multplicity_i
        root_mult_array = [root_mult_array ; root_arr(count,1) i];
        count= count +1;
    end
end

% Print the calculated roots and the corresponding multiplicities.
PrintoutRoots('MY METHOD',root_mult_array);

