function [f_roots,g_roots] = BuildPolynomial(m,intvl_low,intvl_high,seed)
    % ignores the interval size given on input
    % does not consider close roots, which may be shared in f and g, but not
    % considered in d.
    
    a = intvl_low;
    b = intvl_high;
    rng(seed)
    
    
    
    format long
    % Get the multiplicity structure of the roots of the polynomial f
    % t. t = t1 + t2 + ... + t_{r}
    % We want more lower multiplicity roots. so skew this way.
    prob_arr = zeros(1,m);
    for i = 1:1:m
        prob_arr(i) = i./ nchoosek(m+1,2);
    end
    prob_arr = fliplr(prob_arr);
    
    
    
    % Get the multiplicity structure of d.
    total = 0;
    i = 1;
    while total < m
        r = rand;
        prob = prob_arr;
        x = sum(r >= cumsum([0, prob]));
        if (total + x) <= m
            mult_arr_f(i) = x;
            total = total + x;
            i = i+1;
        end
    end
    
    
    
    % get the number of roots of d
    num_roots_f = length(mult_arr_f);
    
    % get a set of unique roots for the polynomials f g and d.
    % the 1000 and 1000 contain the roots to the unit interval
    detail = 100;
    format 'long';
    
    
    roots = a + randperm(detail,num_roots_f)./(detail./(b-a));
    roots_f = roots(1:num_roots_f);
    
    % 
    f_roots = [roots_f' mult_arr_f'];

end