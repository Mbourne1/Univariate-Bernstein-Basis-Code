
m = 5;
n = 7;
k = 2;

A_star = diag(1./(m + n - k : -1 : 1));

B_star = blkdiag(diag(n-k:-1:1), diag(m-k:-1:1));




S = rand(m + n - k + 1, m + n - 2*k + 2);

[Qs, Rs] = qr(S)




