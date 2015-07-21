function u_roots = getQuotient(f_roots,d_roots)
% get the roots of quotient polynomial given the roots of polynomial f, and
% the roots of polynomial d, where d is the GCD of f and g.

num_roots_f_x = size(f_roots,1);
u_roots = [];
for i = 1:1:num_roots_f_x
    % get the root
    root = f_roots(i,1);
    % get multiplicity
    mult_f = f_roots(i,2);
    
    
    % Look up the root in roots of d
    if ~isempty(find(d_roots(:,1) == root));
        % get the row on which we find the root
        [row_d,~] = find(d_roots(:,1) == root);
        
        % get the multiplicity of the root
        mult_d = d_roots(row_d,2);
        
        % subtract multiplicty in d to obtain multiplicity in quotient
        % polynomial u
        mult_u = mult_f - mult_d;
        
        % add the root and its multiplicity to the set of roots for
        % quotient polynomial u
        if mult_u > 0
            u_roots = [u_roots; root mult_u];
        end
        
    else
        
        u_roots = [u_roots; root mult_f];
    end
end