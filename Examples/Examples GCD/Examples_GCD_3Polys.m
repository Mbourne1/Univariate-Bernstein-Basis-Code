function [fx, gx, hx, dx, ux, vx, wx] = Examples_GCD_3Polys(ex_num)
% Inputs. 
%
% ex_num : Example number
%
% Outputs.
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% hx : (Vector) Coefficients of the polynomial h(x)
%
% dx : (Vector) Coefficients of the polynomial d(x), the GCD of f(x) and g(x)
%
% ux : (Vector) Coefficients of the polynomial u(x), given by f(x)/d(x)
%
% vx : (Vector) Coefficients of the polynomial u(x), given by g(x)/d(x)
%
% wx : (Vector) Coefficients of the polynomial u(x), given by h(x)/d(x)



        
[fx, gx, hx, dx, ux, vx, wx] = Examples_GCD_FromCoefficients_3Polys(ex_num);
        
        


end



