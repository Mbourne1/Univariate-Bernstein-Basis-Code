function [vFactors,vMult] = Examples_Deconvolution(ex_num)

% Input f_{i} polynomials
x = sym('x');

% Set example number

switch ex_num
    case '1'
        
        % Create set of factors whose multiplicities are defined in vMult
        vFactors(1,1) = (x+0.017746571505);
        vFactors(2,1) = (x-0.529678501354);
        
        % Set multiplicity of each factor
        vMult = [20 ; 40];
        
    case '2'
        
        % Create set of factors whose multiplicities are defined in vMult
        vFactors(1,1) = (x-0.1259);
        vFactors(2,1) = (x-0.2789);
        vFactors(3,1) = (x-0.389);
        vFactors(4,1) = (x-0.4213);
        vFactors(5,1) = (x-0.5432);
        vFactors(6,1) = (x+0.7923);
        
        % Set multiplicitiy of each factor
        vMult = [ 1; 3; 4; 4; 5; 12 ];
        
    case '3'
        vFactors(1,1) = (x-2);
        vFactors(2,1) = (x-3.2789);
        vFactors(3,1) = (x-1.589);
        vMult = [2; 4; 12 ];
        
    case '4'
        vFactors(1,1) = (x-0.56897);
        vFactors(2,1) = (x+1.24672);
        vFactors(3,1) = (x+0.56921);
        vMult = [3; 6; 9];
        
    case '5'
        
        % Create set of factors whose multiplicities are defined in vMult
        vFactors(1,1) = (x - 0.246512);
        vFactors(2,1) = (x - 1.214654);
        vFactors(3,1) = (x + 0.567890);
        vFactors(4,1) = (x + 0.214654);
        % Set multiplicity of each factor
        vMult = [2; 5; 7; 12];
        
    case '6'
        
        % Create set of factors whose multiplicities are defined in vMult
        vFactors(1,1) = (x - 0.246512);
        vFactors(2,1) = (x - 1.214654);
        vFactors(3,1) = (x + 0.567890);
        vFactors(4,1) = (x + 0.214654);
        % Set multiplicity of each factor
        vMult = [2; 5 ; 7 ; 12];
        
    case '7'
        
        % Create set of factors whose multiplicities are defined in vMult
        vFactors(1,1) = (x-2);
        vFactors(2,1) = (x-3.2789);
        vFactors(3,1) = (x-1.589);
        vFactors(4,1) = (x-0.7213);
        vFactors(5,1) = (x-1.5432);
        vFactors(6,1) = (x+5.72);
        
        % Set multiplicitiy of each factor
        vMult = [ 1; 3; 4 ;4; 5; 12 ];
end