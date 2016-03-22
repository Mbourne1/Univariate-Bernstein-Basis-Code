function [HCG, H1C1G, H2C2G] = BuildHCG(uw,vw,m,n,t)
% Build the matrix HCG, such that H*C(u,v)*G * d = [f;g]

global APF_BUILD_METHOD

% Build H_{t}C(f,g)G_{t}
switch APF_BUILD_METHOD
    case 'Standard'
        
        H1 = BuildH1(m);
        H2 = BuildH1(n);
        C1 = BuildC1(uw,t);
        C2 = BuildC1(vw,t);
        G = BuildG(t);
        
        HCG = blkdiag(H1,H2)*[C1;C2]*G;
        H1C1G = H1*C1*G;
        H2C2G = H2*C2*G;
        
    case 'Rearranged'
        H1C1G = BuildH1C1G(uw,t);
        H2C2G = BuildH1C1G(vw,t);
        HCG = [H1C1G ; H2C2G ];
        
    otherwise
        error('error : Build method is either (Standard) or (Rearranged)');
end


end