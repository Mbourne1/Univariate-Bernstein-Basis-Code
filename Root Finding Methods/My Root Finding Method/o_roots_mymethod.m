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

global SETTINGS

% Initialise an iteration counter for looping through GCD computations.
ite = 1;

% Initialise an array 'q' which stores the gcd outputs from each gcd
% calculation
arr_fx{1,1} = fx;

% let degree_vector store the degrees corresponding to the array of
% GCDs stored in q.
vDeg_arr_fx(1) = GetDegree(arr_fx{1});

% Let theta_vec store all theta values used in each iteration.
vTheta(1) = 1;

% Get the number of distinct roots of f_{1}. Since this is unknown at this
% time, set number of distinct roots to be m_{1} = deg(f_{1}).
vNumber_Distinct_Roots(1) = GetDegree(arr_fx{1});


rank_range = [-8 0];


% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDeg_arr_fx(ite) > 0
    
   
    % if degree of f_{i} is greater than one
    if vDeg_arr_fx(ite) > 1
        
        
        fprintf([mfilename ' : ' sprintf('Compute f_{%i} the GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite-1,ite-1)]);
        
        % Get upper and lower bounds of the GCD Computation.
        % M_{i+1} > M_{i} - d_{i-1}
        try
            
            lowerLimit_t = vDeg_arr_fx(ite) - vNumber_Distinct_Roots(ite-1);
            if lowerLimit_t < 0
                lowerLimit_t = 0;
            end
            upperLimit_t = vDeg_arr_fx(ite) - 1;
            
            fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite, lowerLimit_t)]);
            fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n\n', ite, upperLimit_t)]);
            
        catch
            
            lowerLimit_t = 1;
            upperLimit_t = vDeg_arr_fx(ite)-1;
            
        end
    
        
        
        
        % Get GCD of f(x) and f'(x)
        [arr_fx{ite,1}, ~, arr_fx{ite+1,1}, ~,~,~,~, vDeg_arr_fx(ite+1), rank_range] = ...
            o_gcd_2Polys_mymethod( arr_fx{ite}, Bernstein_Differentiate(arr_fx{ite}), [lowerLimit_t, upperLimit_t], rank_range);
        
        
        % Get number of distinct roots of f(ite)
        vNumber_Distinct_Roots(ite) = vDeg_arr_fx(ite) - vDeg_arr_fx(ite+1);
        
        fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite, vDeg_arr_fx(ite+1))]);
        fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,vNumber_Distinct_Roots(ite))]);

        % increment iteration number.
        ite = ite+1;
        
        
    elseif vDeg_arr_fx(ite) == 1
        
        % if m=1, then n = 0, GCD has maximum degree 0.
        fprintf([mfilename ' : ' 'Only one subresultant exists \n'])
        dx = 1;
        
        
        vDeg_arr_fx(ite+1) = length(dx)-1;
        arr_fx{ite+1} = dx;
        arr_ux{ite} = arr_fx{ite};
        ite = ite+1;
        
        break;
        
        
    end
end


% Get the degree structure of the polynomials h_{i}
vDeg_arr_hx = diff([vDeg_arr_fx]);

% Get the degree structure of the polynomials w_{i}
vDeg_arr_wx = diff([vDeg_arr_hx 0]);


% #########################################################################
% This section performs deconvolution by structured matrix method.


for i = 1:1:length(arr_fx)
    arr_fx{i} = arr_fx{i}./ arr_fx{i}(1);
end


% Deconvolve the set of polynomials f_{i}(x) to get the set of polynomials
% h_{i}(x)

switch SETTINGS.GET_HX_METHOD
    case 'From Deconvolution'
        
        fprintf([mfilename ' : ' sprintf('Deconvolving f_{i}(x) : %s \n',SETTINGS.DECONVOLUTION_METHOD_HX)])
        arr_hx = Deconvolve_Set(arr_fx, SETTINGS.DECONVOLUTION_METHOD_HX);
        
    case 'From ux'
        
        arr_hx = arr_ux;
        
    otherwise
        
        error('err');
        
end

% Get number of polynomials in the array h_{i}(x)
[nPolys_arr_hx] = size(arr_hx,1);



if nPolys_arr_hx == 1
    
    % if number of cols in h1 is only 1, then do not perform second set of
    % deconvolutions, since only one entry in h1.
    
    % Initialise an empty root array
    arr_roots = [];
    
    % Get the polynomial a(w)
    wx = arr_hx{1};
    
    % Normalise the coefficients of a(w)
    wx = wx./wx(1);
    
    % Get the degree of a(w)
    deg_ax = GetDegree(arr_hx{1});
    
    % get aw including its binomial coefficients.
    ax_bi = GetWithBinomials(wx);
    
    % get the roots in terms of z^{i} where z^{i} =(\frac{y}{(1-y)})^{i}.
    % Note we must flip the coefficients to conform with input format of
    % MATLAB roots function.
    rt_wrt_z = roots(flipud(ax_bi));
    
    % Get the roots in terms of y.
    roots_wrt_x = [rt_wrt_z./(1.+rt_wrt_z) ]; %Edit 27/07
    
    % Initialise a vector of ones
    one = ones(length(roots_wrt_x),1);
    
    % Get the roots with respect to y, and their multiplicities all set
    % to one.
    roots_wrt_x = [roots_wrt_x one];
        
    % Add the roots to the array of roots
    arr_roots = [arr_roots ; roots_wrt_x];
else
    
  
    % Deconvolve the second set of polynomials h_{i}(x) to obtain the set of
    % polynomials w_{i}(x)
    arr_wx = Deconvolve_Set(arr_hx, SETTINGS.DECONVOLUTION_METHOD_WX);
    
    
    % w1 yields the simple, double, triple roots of input polynomial f.
    % w1{i} yields the roots of multiplicity i.
    
    % set the w1{max} = h1{max}
    arr_wx{ite-1,1} = arr_hx{ite-1, 1};
    
    % Get number of entries in w1
    [nEntries_arr_wx] = size(arr_wx, 1);
    
    % initialise an empty set
    arr_roots = [];
    
    % for each polynomial w_{i}(x)
    for i = 1:1:nEntries_arr_wx
        
        
        % if the polynomial of said multiplicity is of length 2, degree 1, then
        % only one root exists for this multiplicity. Add it to the list wp1.
        if (length(arr_wx{i}) == 2)

            % Get the polynomial a(w), whose simple roots have multiplicty 
            % i in the polynomial f which we started with.
            wx = arr_wx{i};
            
            % Normalise the polynomial coefficients
            wx = wx./wx(1);
            
            % Convert to power form, so that coefficients are in terms of y^{i}
            % rather than (1-y)^{m-i}y^{i}.
            wx_pwr = [wx(1,:) ; wx(2,:)-wx(1,:)];
            
            % Obtain the root in terms of y, and set multiplicity to one.
            a_root = [-wx_pwr(1,:)./wx_pwr(2,:) wx_pwr(2,:)./wx_pwr(2,:)];
            
            % Add the root to the [root, mult] matrix
            arr_roots = [arr_roots ; a_root];
            
        elseif (length(arr_wx{i}) > 2)
            % The given multiplicity contains more than one root, such that
            % number of coefficients in greater than 2, use MATLAB roots
            % function to find roots.
            % display('Multiplicity contains more than one root. ')
            
            % get the polynomial a(w), whose roots have multiplicity i, in bernstein
            % form.
            wx = arr_wx{i};
            
            % Normalise the polynomial coefficients
            wx = wx./wx(1);
                       
            % get aw including its binomial coefficients.
            ax_bi = GetWithBinomials(wx);
            
            % get the roots in terms of z^{i} where z^{i} =(\frac{y}{(1-y)})^{i}.
            rt_wrt_z = roots(flipud(ax_bi));
            
            % get the root in terms of y
            roots_wrt_x = [rt_wrt_z./(1.+rt_wrt_z) ]; %Edit 27/07
            
            % Initialise a vector of ones
            one = ones(length(roots_wrt_x),1);
            
            % get the roots with respect to y, and their multiplicities all set
            % to one.
            roots_wrt_x = [roots_wrt_x one];
            
            % add the roots to the array of roots
            arr_roots = [arr_roots ; roots_wrt_x];
        end
    end
end


% Obtaining multiplicities of the calculated roots
% create a matrix where the first column contains the multiplicities, and
% the second column contains the number of roots of that multiplicity
nPolys_wi = length(vDeg_arr_wx);
mat = [(1:1:nPolys_wi)' vDeg_arr_wx'];

count = 1;
root_mult_array = [];
for i = 1:1:size(mat,1)
    % Get the number of roots of multiplicity i
    nRoots_of_Multplicity_i = mat(i,2);
    
    for j = 1:1:nRoots_of_Multplicity_i
        root_mult_array = [root_mult_array ; arr_roots(count,1) i];
        count = count +1;
    end
end


end
