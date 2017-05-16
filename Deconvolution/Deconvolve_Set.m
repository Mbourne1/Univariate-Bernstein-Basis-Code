function arr_hx = Deconvolve_Set(arr_fx, DECONVOLVE_METHOD)
% Performs a series of n-1 deconvolutions over a set of n polynomials,
% where each polynomial f_{i}(x) appears in two deconvolutions, and
% h_{i}(x) = f_{i}(x)/f_{i+1}(x)
%
%
% % Inputs
%
% arr_fx : (Array of Vectors) Each cell of the array contains coefficients 
% of the polynomials f_{i}(x) 
%
% % Outputs
%
% arr_hx : (Array of Vectors) Each cell of the array contains coefficients
% of the polynomial h_{i}(x) = f_{i}(x)/f_{i+1}(x)




% Set Deconvolution Method
%   Separate : Use standard non batch method
%   Batch  : Use batch deconvolution
%   Batch With STLN
%   Batch Constrained
%   Batch Constrained With STLN
%


switch DECONVOLVE_METHOD
    case 'Separate'
        
        % Deconvolve independent method
        arr_hx = Deconvolve_Separate(arr_fx);
        
    case 'Batch'
        
        % Deconvolve Batch Method
        arr_hx = Deconvolve_Batch(arr_fx);
        
    case 'Batch With STLN'
        
        arr_hx = Deconvolve_Batch_With_STLN(arr_fx);
        
    case 'Batch Constrained'
        
        % Get number of polynomials in batch
        nPolys_fx = size(arr_fx,1);
        
        % Get the degree of polynomials f_{i}(x)
        vDegt_fx = zeros(nPolys_fx,1);
        for i = 1:1:nPolys_fx
            vDegt_fx(i) = GetDegree(arr_fx{i});  
        end
                
        % Get the degree structure of the polynomials h_{i}
        vDeg_hx = diff(vDegt_fx);
        
        % Get the degree structure of the polynomials w_{i}
        vDeg_wx = diff([vDeg_hx; 0]);
        
        % Get the vector of multiplicity structur of factors of f_{0}(x)
        vMult = find(vDeg_wx~=0);
        
        % Get array of polynomials h_{i}(x)
        arr_hx = Deconvolve_Batch_Constrained(arr_fx,vMult);
        
    case 'Batch Constrained With STLN'
        
        % Get number of polynomials in batch
        nPolys_fx = size(arr_fx,1);
        
        % Get the degree of polynomials f_{i}(x)
        vDegt_fx = zeros(nPolys_fx,1);
        for i = 1:1:nPolys_fx
            vDegt_fx(i) = GetDegree(arr_fx{i});  
        end
                
        % Get the degree structure of the polynomials h_{i}
        vDeg_hx = diff(vDegt_fx);
        
        % Get the degree structure of the polynomials w_{i}
        vDeg_wx = diff([vDeg_hx; 0]);
        
        % Get vector of multiplicity of each factor in f_{0}(x)
        vMult = find(vDeg_wx~=0);
        
        % Compute the polynomials h_{i}(x)
        arr_hx = Deconvolve_Batch_Constrained_With_STLN(arr_fx,vMult);
        
    otherwise
        err_msg = sprintf(...
            [
            'SETTINGS.DECONVOLVE_METHOD must be one of the following:\n'...
            '\t*Separate \n '...
            '\t*Batch \n '...
            '\t*Batch STLN \n '...
            '\t*Batch Constrained\n '...
            '\t*Batch Constrained STLN\n'...
            ]);
        error(err_msg);
end




end
