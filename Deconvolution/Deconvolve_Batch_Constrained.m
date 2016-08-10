function [arr_hx] = Deconvolve_Batch_Constrained(arr_fx,vMult)
% Let vMultiplicities be a vector containing the multiplicities of the
% roots of f_{0}(x). vMult = [m_{1}, m_{2} ,..., m_{n}]
% The division f_{0}/f_{1},...,f_{m_{1}-1} / f_{m_{1}} all have the same solution
%
%
% % Inputs.
%
% arr_fx : Array of polynomials f(x)
%
% vMult : Multiplicities of the factors of f_{0}(x) in ascending order.
%
% % Outputs.
%
% arr_hx : Array of polynomials h(x) where h_{i}(x) = f_{i-1}(x) / f_{i}(x) 


% %
% %
% Build LHS Matrix

LHS_Matrix = BuildDTQ(arr_fx,vMult);


% %
% %
% Build the RHS vector



% RHS vector consists of f_{1},...,f_{m_{i}} where m_{i} is the highest
% degree of any root of f_{0}(x).
RHS_vec = [];
for i = 1:1:length(arr_fx)-1
    RHS_vec = [RHS_vec ; arr_fx{i}];
end

x = SolveAx_b(LHS_Matrix,RHS_vec);

x_temp = x;
%for i = 1:1:length(vMult)
unique_vMult = unique(vMult);

for i = 1:1:length(unique_vMult)
    mult = unique_vMult(i);
    deg = GetDegree(arr_fx{mult}) - GetDegree(arr_fx{mult+1});
       
    arr_px{i} = x_temp(1:deg+1);
    x_temp(1:deg+1) = [];
end


% % 
% %
% Get the polynomials p_{i}(x) repeated to give the set of polynomials
% h_{i}(x).

nEntries_px = length(arr_px);

count = 1;
for i = 1:1:nEntries_px
        
    if i == 1
        nReps = unique_vMult(i);
    else
        nReps = (unique_vMult(i) - unique_vMult(i-1));
    end
        
    for j = 1:1:nReps
        arr_hx{count} = arr_px{i};
        count = count + 1; 
    end
    
end

% %
% %
% Get Residual
residual = RHS_vec - (LHS_Matrix*x);
display(norm(residual));


end

function LHS_Matrix = BuildDTQ(arr_fx,vMult)
% % 
% %
% Build the LHS Coefficient matrix

% For each distinct hx, build the partition of the Matrix C(f_{m_{i}+1},...,f_{m_{i+1}})
% C(f_{1},...f_{m1}), C(f_{m1+1},...,f_{m2}),...
% C(f_{1},...,f_{m1}) = [T(f1) ; T(f2) ; ... T(fm1)]

% Get number of distinct polynomials in h_{i} 
nDistinct_hx = length(vMult);

for i = 1:1:nDistinct_hx
    
    if i > 1 
        old_mult = vMult(i-1);
    else % No previous multiplicity
        old_mult = 0;
    end
    
    new_mult = vMult(i);
    
    Cf{i} = [];
    
    % for each polynomial f_{i} in the interval f_{m_{i-1}+1}...f_{m_{i}}
    for j = (old_mult+1+1) : 1 : (new_mult+1)
        
        % Get the degree of the previous polynomial f_{i-1}(x)
        fx_prev = arr_fx{j-1};
        deg_fx_prev = GetDegree(fx_prev);
        
        % Get the degree of the current polynomial f_{i}(x)
        fx = arr_fx{j};
        deg_fx = GetDegree(fx);
        
        % Get the degree of the polynomial h_{i} = f_{i-1}/f_{i}
        deg_hx = deg_fx_prev - deg_fx;
        
        % Build the Cauchy like matrix T_{m_{i} - m_{i-1}}(f_{i})
        Tf{j} = BuildT1(fx, deg_hx);
        
        D{j} = BuildD(deg_fx,deg_hx);
        
        % Stack beneath all other T_{f} which are multiplied by [_{i}(x)
        Cf{i} = [Cf{i} ; D{j}*Tf{j}];
    end
    
    Q{i} = BuildQ1(deg_hx);
    
    DTQ{i} = Cf{i} * Q{i};
    
end

LHS_Matrix = blkdiag(Cf{:});
LHS_Matrix = blkdiag(DTQ{:});

end