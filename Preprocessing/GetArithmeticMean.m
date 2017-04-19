function lambda = GetArithmeticMean(fx, n_k)
%
% % Inputs
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% n_k : (Int) Degree of polynomial v(x)
%
% % Outputs
%
% lambda : (Float) Arithmetic mean of non-zero entries of the partition
% T_{n-k}(f(x)).

global SETTINGS



switch SETTINGS.SYLVESTER_BUILD_METHOD 
    % Possible Values.
    % T : 
    % DT :
    % DTQ :
    % TQ : 
    % DTQ Rearranged Denom Removed :
    % DTQ Rearranged :
    
    
    case 'DTQ'
    
        m = GetDegree(fx);
        
        lambda = ((m+n_k+1) / ((m+1)*(m+1)*(n_k+1))) * sum(abs(fx));
        
        
        %Tf = BuildSubresultant_Partition_2Polys(fx, n_k);
        %lambda2 = mean(Tf(Tf~=0));
        %display(lambda)
        
    case {'T', 'DT','TQ', 'DTQ Rearranged Denom Removed','DTQ Rearranged'}
        
        Tf = BuildSubresultant_Partition_2Polys(fx, n_k);
        lambda = mean(Tf(Tf~=0));
        
    otherwise 
       
        error('Not a valid Sylvester Matrix build method')
end

end