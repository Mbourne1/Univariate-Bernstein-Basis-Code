function [HCG, H1C1G, H2C2G] = BuildHCG(uw,vw,m,n,t)
% Build the matrix HCG, such that H*C(u,v)*G * d = [f;g]
%
% uw : Matrix of coefficients of the polynomial u(w,w) 
%
% vw : Matrix of coefficients of the polynomial v(w,w)
%
% m : Degree of polynomial f(x,y)
%
% n : Degree of polynomial g(x,y)
%
% t : Degree of GCD d(x,y)


global SETTINGS

% Build H_{t}C(f,g)G_{t}
switch SETTINGS.APF_BUILD_METHOD
    case 'Standard'
        
        H1 = BuildH1(m);
        H2 = BuildH1(n);
        
        C1 = BuildT1(uw,t);
        C2 = BuildT1(vw,t);
        
        G = BuildG(t);
        
        HCG = blkdiag(H1,H2)*[C1;C2]*G;
        H1C1G = H1*C1*G;
        H2C2G = H2*C2*G;
        
    case 'Rearranged'
        H1C1G = BuildH1C1G(uw,t);
        H2C2G = BuildH1C1G(vw,t);
        HCG = [H1C1G ; H2C2G ];
        
    otherwise
        error('error : SETTINGS.APF_BUILD_METHOD method is either (Standard) or (Rearranged)');
end


end