function [] = Test_Deconvolution
% Test the different methods of deconvolving the polynomials f_{i}(x), to
% form the set of polynomials h_{i} where h_{i} = f{i}/f_{i+1}
%
%
%
%
%

% Set settings pertaining to this test

global SETTINGS
SETTINGS.PLOT_GRAPHS = 'y';
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-15;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 100;
SETTINGS.BOOL_LOG = 'n';
SETTINGS.BOOL_DENOM_SYL = 'y';
SETTINGS.SYLVESTER_BUILD_METHOD = 'Standard';

% Input f_{i} polynomials
x = sym('x');

% Set example number 
ex_num = '1';

switch ex_num
    case '1'
        
        % Create set of factors whose multiplicities are defined in vMult
        factor(1) = (x - 0.246512);
        factor(2) = (x - 1.214654);
        factor(3) = (x + 0.567890);
        factor(4) = (x + 0.214654);
        % Set multiplicity of each factor
        vMult = [2, 5 , 7 , 12];
        
    case '2'
                                
        % Create set of factors whose multiplicities are defined in vMult
        factor(1) = (x-2);
        factor(2) = (x-3.2789);
        factor(3) = (x-1.589);
        factor(4) = (x-0.7213);
        factor(5) = (x-1.5432);
        factor(6) = (x+5.72);
        
        % Set multiplicitiy of each factor
        vMult = [ 1 3 4 4 5 12 ];
end

% Get highest power of any factor
highest_pwr = max(vMult);

% %
% %
% Generate polynomials f_{0}(x) ,..., f_{m}(x) = 1. Where each f_{i+1}(x) is
% the f_{i+1} = GCD(f_{i},f'_{i}).
for i = 0:1:highest_pwr
    
    % Get the multiplicities of the roots of f_{i+1}
    mults = ((vMult - i) + abs(vMult-i)) ./2;
    
    % Get the symbolic polynomial f_{i+1}
    arr_sym_f{i+1} = prod(factor.^(mults));
    
    % Get the degree of polynomial f_{i+1}(x)
    vDeg_f(i+1) = double(feval(symengine, 'degree', (arr_sym_f{i+1})));
end


% Get the degree structure of the polynomials h_{i} where h_{i} =
% f_{i-1}(x)/f_{i}(x)
vDeg_arr_hx = diff(vDeg_f);

% Get the degree structure of the polynomials w_{i} where w_{i} =
% h_{i-1}/h_{i}
vDeg_arr_wx = diff([vDeg_arr_hx 0]);

% Get the multiplicities of the roots.
vMultiplicities = find(vDeg_arr_wx~=0);

% Get the sequence of polynomials h_{i}(x) in symbolic form
sym_arr_h = cell(length(arr_sym_f)-1,1)
for i = 1:1:length(arr_sym_f)-1
    sym_arr_h{i} = arr_sym_f{i} / arr_sym_f{i+1};
end

% %
% %
% Get coefficients vectors of f_{i}(x) and h_{i}(x)
nPolys_arr_fx = length(arr_sym_f);
nPolys_arr_hx = length(arr_sym_f) - 1;

arr_fx = cell(nPolys_arr_fx,1);
arr_hx = cell(nPolys_arr_hx,1);

for i = 1:1:nPolys_arr_fx
    if i <= nPolys_arr_hx
        arr_fx{i,1} = sym2poly(arr_sym_f{i})';
        arr_hx{i,1} = sym2poly(sym_arr_h{i})';
    else
        arr_fx{i,1} = 1;
    end
    
end

% %
% %
% Convert to Bernstein Basis
for i = 1:1:length(arr_fx)
    arr_fx{i,1} = PowerToBernstein(arr_fx{i,1})
end

for i = 1:1:length(arr_hx)
    arr_hx{i,1} = PowerToBernstein(arr_hx{i,1})
end



% %
% %
% %
% Add noise to the coefficients of f_{i}(x)

arr_fx_noisy = cell(nPolys_arr_fx,1);

% Set upper and lower noise levels
emin = 1e-8;
emax = 1e-12;
for i = 1:1:nPolys_arr_fx
    arr_fx_noisy{i,1} = Noise(arr_fx{i},emin,emax);
end




% -------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution batch method with constraints

LineBreakLarge()
fprintf('Deconvolve with the constraint that all h_{i} for i between m_{j} and m_{j+1} are equal \n')
fprintf('Staggered Staircase Method \n\n')

arr_hx_test_1 = Deconvolve_Batch_Constrained(arr_fx_noisy,vMultiplicities);

for i = 1:1:length(arr_hx_test_1)
    %display(arr_hx_test_1{i} );
end

% Compare each computed h{i} with actual h_{i}
for i = 1:1:length(arr_hx_test_1)
    
    exact = arr_hx{i};
    comp = arr_hx_test_1{i};
    
    err_measure = norm(exact - comp) ./ norm(exact);
    fprintf([mfilename ': ' sprintf('%i : Error : %2.4e \n', i,err_measure)]);
end

% -------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution batch method with constraints and low rank approx

% LineBreakLarge()
% fprintf('Deconvolve with the constraint that all h_{i} for i between m_{j} and m_{j+1} are equal \n')
% fprintf('With STLN \n')
% fprintf('Staggered Staircase Method \n\n')
% 
% arr_hx_test_3 = Deconvolve_Batch_Constrained_With_STLN(arr_fx_noisy,vMultiplicities);
% 
% for i = 1:1:length(arr_hx_test_3)
%     %display(arr_hx_test_1{i} );
% end
% 
% % Compare each computed h{i} with actual h_{i}
% for i = 1:1:length(arr_hx_test_3)
%     
%     exact = arr_hx{i};
%     comp = arr_hx_test_3{i};
%     
%     err_measure = norm(exact - comp) ./ norm(exact);
%     fprintf([mfilename ': ' sprintf('%i : Error : %2.4e \n', i,err_measure)]);
% end

%--------------------------------------------------------------------------
% %
% %
% %
% Testing standard deconvolution batch method
LineBreakLarge();
fprintf('Deconvolution Batch (Staircase Method)\n')

arr_hx_test_3 = Deconvolve_Batch(arr_fx_noisy);
for i = 1:1:length(arr_hx_test_3)
    %display(arr_hx_test_2{i});
end

% Compare each computed h{i} with actual h_{i}
for i = 1:1:length(arr_hx_test_3)
    
    exact = arr_hx{i};
    comp = arr_hx_test_3{i};
    
    err_measure = norm(exact - comp) ./ norm(exact);
    fprintf([mfilename ': ' sprintf('%i : Error : %2.4e \n', i,err_measure)]);
end



%--------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution
LineBreakLarge();
fprintf('Deconvolve - Each deconvolution is independent \n');
fprintf('Deconvolution Separate \n');


arr_hx_test_4 = Deconvolve_Separate(arr_fx_noisy);

for i = 1:1:length(arr_hx_test_4)
    %display(arr_hx_test_3{i});
end
% Compare each computed h{i} with actual h_{i}
for i = 1:1:length(arr_hx_test_4)
    
    exact = arr_hx{i};
    comp = arr_hx_test_4{i};
    
    err_measure = norm(exact - comp) ./ norm(exact);
    fprintf([mfilename ': ' sprintf('%i : Error : %2.4e \n', i,err_measure)]);
end


end