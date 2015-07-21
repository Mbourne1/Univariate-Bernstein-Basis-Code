function [f_noisy]=Noise(f,el,eu)

% Add noise to the coefficients of polynomial f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs

% f :

% el : signal to noise low limit

% eu : signal to noise upper limit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the degree of input polynomial f
m = length(f) - 1;

switch nargin
    case 2 % Only one noise is specified, set upper = lower
        
        
        rp = (2*rand(1,m+1))-ones(1,m+1);
        s = rp*el;
        
        noisevector = f.*s;
        f_noisy = f + noisevector;
        
        
    case 3 % Specified upper and lower bound of noise
        
        y = (2*rand(m+1,1))-ones(m+1,1);
        s = eu *ones(m+1,1) -  y.*(eu-el);
        noisevector = f'.*s;
        f_noisy = f + noisevector';
end

