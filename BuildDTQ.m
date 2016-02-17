function [DTQ] = BuildDTQ(fx_n,gx_n,alpha,theta,t)

global bool_sylvesterBuildMethod

[r,~] = size(fx_n);
m = r - 1;

[r,~] = size(gx_n);
n = r - 1;


switch bool_sylvesterBuildMethod
    case 'standard'
        % Build Matrices D,T and Q.
        D = BuildD(m,n,t);
        T = BuildT(fx_n,gx_n,alpha,theta,t);
        Q = BuildQ(m,n,t);
        DTQ = D*T*Q;
        
    case 'rearranged'
        DTQ = BuildDTQ_ElementWise(fx_n,gx_n,alpha,theta,t);
    otherwise
        error('Error : Build Method is either (standard) or (rearranged)')
end

end