function Y = BuildDYNew(m,n,t,idxMinCol,x_ls_wrt_w,alpha,theta)
%
% Inputs.
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% t : Degree of polynomial d(x)
%
% idxMinCol : Index of minimum column
%
% x_ls : 
%
% alpha : 
%
% theta :

xa = x_ls_wrt_w(1:idxMinCol-1) ;
xb = x_ls_wrt_w(idxMinCol:end) ;
x_ls_wrt_w = [xa; 0 ;xb] ;% Insert zero into vector

% Get number of cols in left partition
c_leftmatrix = (n-t+1);

% Get number of cols in right partition
c_rightmatrix = (m-t+1);


X1 = x_ls_wrt_w(1:c_leftmatrix);

X2 = x_ls_wrt_w(c_leftmatrix+1:c_leftmatrix+c_rightmatrix);


%Build Empty Sylvester Matrix
% First half
Y1 = zeros(m+n-t+1,m+1);

% Second half
Y2 = zeros(m+n-t+1,n+1);

% Build the first half of the matrix DY
% For each column j = 0,...,
for j=0:1:m
    % For each row i = j,...,
    for i=j:1:j + length(X1)-1
        Y1(i+1,j+1) = ...
            X1(i-j+1) .*(theta^(j)) ...
            .* nchoosek(m,j) ...
            .* nchoosek(n-t,i-j)./ ...
            nchoosek(m+n-t,i);
    end
end



% Build the second half of the matrix DY
% for each column j = 0,...,n
for j=0:1:n
    % for each row i = j,...,
    for i=j:1:j + length(X2)-1
        Y2(i+1,j+1)= ...
            X2(i-j+1).*(theta^(j)) .* ...
            nchoosek(n,j) .*...
            nchoosek(m-t,i-j)./...
            nchoosek(m+n-t,i);
    end

end



Y=[Y1,alpha.*Y2];

end