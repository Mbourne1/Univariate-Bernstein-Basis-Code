function C = Cauchy(c,t) 
%% Create the Cauchy Matrix for Bernstein polynomials.
%Let c be the quotient polynomial whose cauchy matrix is to be built
%Let t be the degree of the divisor
    
    x=length(c)-1;
    m = x+t; % The degree of either f or g
    
    %Let c be a cauchy matrix of the coefficients

    %for each column
    for k=1:1:t+1
        %for each row
        for h=k:1:k+x
            C(h,k)= c(h-k+1);
        end
    end

    
end
