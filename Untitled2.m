function [] = Untitled2

x = sym('x');

p = 1;
q = 1;

for j = 1:1:10
    
xj(j) = (-1)^j* (j/2);
end
nRoots = length(xj);

for j = 1:1:nRoots
   p = p .* (x - xj(j));
end

for j = 1:1:nRoots
   q = q .* (x-xj(j) +10^(-j)); 
   %q = q .* (x-xj(j)); 
end

display(p)
display(q)

f = fliplr(sym2poly(p))';
g = fliplr(sym2poly(q))';

f = PowerToBernstein(f);
g = PowerToBernstein(g);

display(f)
display(g)

m = GetDegree(f);
n = GetDegree(g);



problemType = 'From';
ex_num = '1';
emin = 0;
emax = 0;
%mean_method = 'Geometric Mean Matlab Method';
%bool_alpha_theta = 'y';
mean_method = 'None';
bool_alpha_theta = 'n';
low_rank_approx_method = 'None';
apf_method = 'None';

SetGlobalVariables(problemType,ex_num,emin,emax,mean_method,...
    bool_alpha_theta,low_rank_approx_method, apf_method);


[f,g,d,u,v] = o_gcd_mymethod(f,g,[0,min(m,n)],'n')

end