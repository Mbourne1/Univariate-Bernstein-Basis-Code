

m = 5;
n = 10;
k = 2;

fx = ones(m+1,1);
gx = ones(n+1,1);

Sk = BuildDTQ(fx,gx,k);

vec = m + n - k : -1 : 1;
A = [diag(1./vec) zeros(m + n - k, 1)];
[Qa,Ra] = qr(A);

display(Qa)
display(Ra)


[Qs, Rs] = qr(Sk);

Sk2 = Sk;
Sk2(end,:) = [];
[Qs2, Rs2] = qr(Sk2)

display(Qs)
display(Rs)

display(Qs2)
display(Rs2)

vec = n - k : -1 : 1;
B1 = [diag(vec) ; zeros(1, n - k )];

vec = m - k : -1: 1;
B2 = [diag(vec) ; zeros(1, m - k )];

B = blkdiag(B1,B2);

[Qb,Rb] = qr(B);
display(Qb)
display(Rb)


fprintf('Experiments \n')

