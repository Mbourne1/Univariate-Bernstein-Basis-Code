function s = Subresultant(m,n)


% start value of k
k = 0;

% Build C1

    C_0 = BuildCauchy(m,n,k)

% Build C1 where k = 1

    
    C_1 = BuildCauchy(m,n,k+1)


    

    C_1b = BuildSubresultant(m,n,k,C_0)

    
   
    

end



        

function C1 = BuildSubresultant(m,n,k,C_0)
% Where C_0 is the previous Cauchy Matrix,
    
    knew = k+1;
% Build Matrix A

    A = [zeros((m+n-2*(knew)+2),1) diag(1./(1:1:(m+n-2*(k)+2))) ]
    
% Build Matrix B


B = [...
        zeros(1,n-k);...
        diag(1:1:(n-k))
        ];
    
% Build C1

    C1 = A * C_0 * B;

    rat = nchoosek(m+n-(k),n-(k)) ./ nchoosek(m+n-(knew),n-(knew));
    
    C1 = C1 .* rat;





end










