function [] = Bisection_method(example_number,emin)
    
    addpath 'Examples'

    [f_root_mult_arr_exact,~] = Root_Examples(example_number);
    
    % Get Coefficients of f(x) in scaled Bernstein form
    f_exact_bi = BuildPolyFromRoots(f_root_mult_arr_exact);
    
    % Get the degree of f(x)
    m = GetDegree(f_exact_bi)
    
    % Get f(x) without binomials, in standard Bernstein form
    f_exact = GetWithoutBinomials(f_exact_bi);

    
    % Set signal noise ratio, note only min is specified, so we have
    % constant signal/noise ratio
    fx = Noise(f_exact,emin);
    
    index = 0;
    for i = 0:0.05:1
        index = index + 1; 
        yy(index) = BernsteinEval(fx,i);
    end
    xx = 0:0.05:1;
    
    
%% Control polygon
% Let Pk be the set of control points, P_0,P_1,...,P_m 
    Pk = [];
   
% Let yk be the set of coefficients
    yk = fx;
    
    a = 0;
    b = 1;

% Pk is the set of control points of the curve
    for i = 0:1:m
        Pk = [Pk ; a + (i/m).*(b-a)     yk(i+1)];
    end
    
    fprintf('The set of control Points')
    Pk
    
    c = 0.21
    Pk_array{1} = Pk
   
    
% Pk_array{i} contains the set of control points P^{i}. 
% subdivision at point c to obtain P_left and P_Right
    
    % obtain each set of control points
    for i = 2:1:m+1
       Pk_array{i} = deCasteljau(c,Pk_array{i-1})  ;
    end
    
    for i = 1:1:length(Pk_array)
       plot(Pk_array{i}(:,1),Pk_array{i}(:,2),'black-'); 
    end
    
    
% From the generated sets of control points obtain set of control points 
% for sub-curve P[a,c]
    new_left_pk = zeros(m+1,2);
    for i = 1:1:m+1
        new_left_pk(i,:) = Pk_array{i}(1,:);
    end
    
% obtain set of control points for sub-curve P[c,b]
    new_right_pk = zeros(m+1,2);
    for i=1:1:m+1
        new_right_pk(i,:) = Pk_array{i}(end,:);
    end
    
% Graph the calculated control points.
    figure(1)
    plot(xx,yy);
    hold on
    x  = Pk(:,1);
    y = Pk(:,2)
    plot(x,y,'linewidth',2);
    
    k = convhull(Pk)
    X = Pk
    
         K = convhulln(X);
         P = [];
         for k = 1:size(K,1)
          Q = X(K(k,:),:).';
          Q1 = Q(:,1);
          R = [P1-P2,Q(:,2:n)-repmat(Q1,1,n-1)];
          s = R\(P1-Q1);
          t = s(1); s(1) = 1-sum(s(2:n));
          if all(s>=0)

           P = [P,(1-t)*P1+t*P2];
          end
         end
    
    
    plot(x(k),y(k),'r-');
    
    
    %plot(Pk(:,1),Pk(:,2),'black-');


    %plot(new_left_pk(:,1),new_left_pk(:,2),'red-s');
    %plot(new_right_pk(:,1),new_right_pk(:,2),'red-o');
    hold off
    
        
    
    
    
    
end


function [xx] =  BernsteinEval(f,c)
%% Evaluate function f at point c \in t
 y = c;
 m = length(f)-1;
 
    for i = 0:length(f)-1
        
        x(i+1) = f(i+1) .* nchoosek(m,i) .* (1-y)^(m-i) .* y^i ;
    end
    
    xx = sum(x);
    
end

function Pk_1 = deCasteljau(c, Pk_0)
% Obtain the next set of control points, given a set of control points Pk_0
% in the deCasteljau subdivision process
% a = start point on horizontal axis
% b = end point on horizontal axis
% c = subdivision point on horizontal axis
% Pk_0 = original control points

% get first control poitn
    a = Pk_0(1,1);
% get last control point
    b = Pk_0(end,1);

    
    T = (c-a)./(b-a);
    
    
% Build set of control points
    Pk_1 = [];
    % for 
    for i = 1:1:length(Pk_0)-1
        %Pk_1(i)
        aa =  (1-T).*Pk_0(i,:) + T.*Pk_0(i+1,:) ;
        Pk_1 = [Pk_1 ; aa];
    
    end
    
    

end
