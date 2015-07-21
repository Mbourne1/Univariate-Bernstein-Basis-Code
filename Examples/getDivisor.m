
function d_roots = getDivisor(f_roots,g_roots)
% get the common roots of polynomials f and g.

% get the number of roots in polynomial f
num_roots_f = size(f_roots,1);

d_roots = [];

% for each root in f_x, check to see if it exists in g_x
for i = 1:1:num_roots_f
    root = f_roots(i,1);
    mult_root_in_f = f_roots(i,2);
    
    % look up the root in g
    if ~isempty(find(g_roots(:,1) == root));
        [row_d,~] = find(g_roots(:,1) == root);
        mult_root_in_g = g_roots(row_d,2);
        mult_root_in_d = min(mult_root_in_f,mult_root_in_g);
        d_roots = [d_roots ; root mult_root_in_d]; 
    end
end


end

