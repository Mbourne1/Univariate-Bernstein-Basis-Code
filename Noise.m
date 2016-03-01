function [f_noisy]=Noise(f,el,eu)
% Add noise to the coefficients of polynomial f(x)
%
%                           Inputs
%
% f  : Column vector of coefficients of polynomial f(x).
%
% el : signal to noise low limit.
%
% eu : signal to noise upper limit
%


% get the degree of input polynomial f
[r,c] = size(f);
m = r -1;

switch nargin
    case 2 % Only one noise is specified, set upper = lower
        
        
        rp = (2*rand(m+1,1))-ones(m+1,1);
        s = rp*el;
        
        noisevector = f.*s;
        f_noisy = f + noisevector;
        
        
    case 3 % Specified upper and lower bound of noise
        
        y = (2*rand(m+1,1))-ones(m+1,1);
        s = eu *ones(m+1,1) -  y.*(eu-el);
        noisevector = f'.*s;
        f_noisy = f + noisevector';
end

