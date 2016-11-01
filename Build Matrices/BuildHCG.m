function [HCG, H1C1G, H2C2G] = BuildHCG(uw,vw,t)
% Build the matrix HCG, such that H*C(u,v)*G * d = [f;g]
%
% uw : Matrix of coefficients of the polynomial u(w,w) 
%
% vw : Matrix of coefficients of the polynomial v(w,w)
%
% t : Degree of GCD d(x,y)


% Build the first part
H1C1G = BuildH1C1G(uw,t);

H2C2G = BuildH1C1G(vw,t);

HCG = ...
    [
        H1C1G;
        H2C2G
    ];

end