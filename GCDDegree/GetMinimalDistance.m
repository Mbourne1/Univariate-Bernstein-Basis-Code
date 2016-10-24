function MinimumDistance = GetMinimalDistanceQR(Sk,method)
% Given a sylvester matrix S_{k}, get the minimum residual 
% Ax = b, x = pinv(A)b => r = b - Ax

% Get number of columns in subresultant S_{k}
[~,nCols_Sk] = size(Sk);


% Initialise a vector to store residuals
vResidual = zeros(1,nCols_Sk);

% For each column in S_{k}(f,g)
for i = 1:1:nCols_Sk
    [Ak,ck]   = RemoveSubresultantColumn(Sk,i);
    
    % Pick method for obtaining residuals
    
    
    switch method
        case 'QR' % Get Residual by QR method
            vResidual(i) = CalculateResidualQR(ck,Ak);
        case 'SVD' % Get Residual by SVD method
            vResidual(i) = CalculateResidualSVD(ck,Ak);
    end
end

MinimumDistance = min(vResidual);

end