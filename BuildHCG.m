function [HCG, H1C1G, H2C2G] = BuildHCG(uw,vw,m,n,t)

global Bool_APFBuildMethod

% Build H_{t}C(f,g)G_{t}
switch Bool_APFBuildMethod
    case 'standard'
        H1 = BuildH1(m);
        H2 = BuildH1(n);
        C1 = BuildC1(uw,t);
        C2 = BuildC1(vw,t);
        G = BuildG(t);
        
        HCG = blkdiag(H1,H2)*[C1;C2]*G;
        H1C1G = H1*C1*G;
        H2C2G = H2*C2*G;
        
    case 'rearranged'
        H1C1G = BuildH1C1G(uw,t);
        H2C2G = BuildH1C1G(vw,t);
        HCG = [H1C1G ; H2C2G ];
        
    otherwise
        error('error : Build method is either (standard) or (rearranged)');
end


end