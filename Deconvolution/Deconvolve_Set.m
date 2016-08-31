function arr_hx = Deconvolve_Set(arr_fx,DECONVOLVE_METHOD)
% Performs a series of d deconvolutions over a set of polynomials,
% where each polynomial g_{i} appears in two deconvolutions.
%
%
% % Inputs
%
% set_f :   Set of input polynomials g(y) to be deconvolved. Each g_{i} has a
%           different number of elements, so set_g is a cell array.
%
% % Outputs
%
% h_{i} = g_{i-1}/g_{i}
%
%



% Set Deconvolution Method
%   single := Use standard non batch method
%   batch  := Use batch deconvolution

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
        
        
        vMult = find(vDeg_wx~=0);
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
        
        
        vMult = find(vDeg_wx~=0);
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
