ex_num = '1';
[fx, gx, dx, ux, vx] = Examples_GCD(ex_num);


m = GetDegree(fx);
n = GetDegree(gx);
k = 2;

A_star = diag(1./(m + n - k : -1 : 1));
I_A = [eye(m + n - k) zeros(m + n -k,1)];

B_star = blkdiag(diag(n-k:-1:1), diag(m-k:-1:1));
I_B = [...
    eye(n - k , n - k ) zeros(n - k , m - k ) ;...
    zeros(1, n - k ) zeros(1, m - k ) ;...
    zeros(m - k , n - k ) eye(m - k , m - k ) ; ...
    zeros(1, n - k ) , zeros(1, m - k )
    ];

S = BuildSubresultant_2Polys(fx,gx,k);

S = I_A * S * I_B;

[Qs, Rs] = qr(S)

[Qs2, Rs2] = qr(A_star * S)


Qs./Qs2

Rs./ Rs2
display('end')




