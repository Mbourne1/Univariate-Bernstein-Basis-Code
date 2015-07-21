function [f_roots,g_roots] = BuildRandomPolynomials2(m,n,t,intvl_low,intvl_high,seed)
    % ignores the interval size given on input
    % does not consider close roots, which may be shared in f and g, but not
    % considered in d.
    
    a = intvl_low;
    b = intvl_high;
    
    format long
    % Get the multiplicity structure of the roots of the GCD d of degree
    % t. t = t1 + t2 + ... + t_{r}
    % We want more lower multiplicity roots. so skew this way.
    prob_arr = zeros(1,t);
    for i = 1:1:t
        prob_arr(i) = i./ nchoosek(t+1,2);
    end
    prob_arr = fliplr(prob_arr);
    
    % Get the multiplicity structure of d.
    total = 0;
    i = 1;
    while total < t
        r = rand;
        prob = prob_arr;
        x = sum(r >= cumsum([0, prob]));
        if (total + x) <= t
            mult_arr_d(i) = x;
            total = total + x;
            i = i+1;
        end
    end
    
    % get the number of roots of d
    num_roots_t = length(mult_arr_d);
    
    
    
    
    % get multiplicity structure for the remaining roots of f
    % initialise a probability vector so that lower multiplicities are
    % preferred.
    prob_arr = zeros(1,m-t);
    for i = 1:1:(m-t)
        prob_arr(i) = i./ nchoosek((m-t)+1,2);
    end
    
    prob_arr = fliplr(prob_arr);
    sum(prob_arr)
    % get the multiplicity structure of the remaining roots of f.
    total = 0;
    i = 1;
    while total < m-t
        r = rand;
        prob = prob_arr;
        x = sum(r >= cumsum([0, prob]));
        if (total + x) <= m-t
            mult_arr_f(i) = x;
            total = total + x;
            i = i+1;
        end
    end
    
    num_roots_f = length(mult_arr_f);
    mult_arr_f = [mult_arr_d mult_arr_f];
    
    
    % initialise a probability vector so that lower multiplicities are
    % preferred
    prob_arr = zeros(1,n-t);
    for i = 1:1:(n-t)
        prob_arr(i) = i./ nchoosek((n-t)+1,2);
    end
    format short
    prob_arr = fliplr(prob_arr);
    
    sum(prob_arr)
    total = 0;
    i = 1;
    while total < n-t
        r = rand;
        prob = prob_arr;
        x = sum(r >= cumsum([0, prob]));
        if (total + x) <= n-t
            mult_arr_g(i) = x;
            total = total + x;
            i = i+1;
        end
    end
    
    num_roots_g = length(mult_arr_g);
    mult_arr_g = [mult_arr_d mult_arr_g];
    
    
    % get a set of unique roots for the polynomials f g and d.
    % the 1000 and 1000 contain the roots to the unit interval
    detail = 100;
    format 'long';
    rng(seed)
    roots = a + randperm(detail,num_roots_t+num_roots_g + num_roots_f)./(detail./(b-a));
    
    
    roots_d = roots(1:num_roots_t);
    roots(1:num_roots_t) = [];
    roots_f = roots(1:num_roots_f);
    roots(1:num_roots_f) = [];
    roots_g = roots(1:num_roots_g);

    
    f_roots = [[roots_d'; roots_f'] mult_arr_f'];
    g_roots = [[roots_d'; roots_g'] mult_arr_g'];
end