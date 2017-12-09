

m = 7;
n = 5;
k = 2;

fx = ones(m + 1,1);
gx = ones(n + 1,1);

D = BuildD_2Polys(m,n-k);
T = BuildT_2Polys(fx,gx,k);

Sk = D * T;

D2 = BuildD_2Polys(m, n-k-1);
T2 = BuildT_2Polys(fx, gx, k+1);
SkPlus1 = D2 * T2;


vec = m + n - k : -1 : 1;
A = [diag(1./vec) zeros(m + n - k, 1)];

vec = ones(n - k,1);
B1 = [diag(vec) ; zeros(1, n - k )];

vec = ones(m - k,1);
B2 = [diag(vec) ; zeros(1, m - k )];

B = blkdiag(B1,B2); 

SkPlus1_test = (m + n - k) .* A * Sk * B

display(Sk)
display(SkPlus1)
display(SkPlus1_test)


fprintf('These two things should be the same \n')
%[Qa,Ra] = qr(A);
[Qs,Rs] = qr((m + n - k) .* A *Sk*B);
[Qb,Rb] = qr(B);
fprintf('Second thing\n')



[Q,R] = qr((m + n - k) .* A * Sk)
[qupdate, rupdate] = qrdelete(Q,R,m + n -(2*k)+2,'col');
[qupdate, rupdate] = qrdelete(qupdate,rupdate,n-k+1,'col');

fprintf('end')

% Check equivalence 
[Qa, Ra] = qr(A);
[Qs, Rs] = qr(Sk);

[Q1, R1] = qr(Qa * Ra * Qs * Rs)

[Q2, R2] = qr(A*Sk)

fprintf('end')


[Q3, R3] = qr(A*Qs)

display(Q3)
display(R3 * Rs)


test1 = A * Sk
test2 = Q3 * R3 * Rs

fprintf('end')


