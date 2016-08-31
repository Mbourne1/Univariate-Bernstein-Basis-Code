
function fx_exact = Examples_Roots_FromCoefficients(ex_num)
switch ex_num
    case '1'
        fx_exact = ...
            [
            -0.9865
            2.2398
            2.8950
            1.9092
            -0.1477
            ];
    otherwise
        error('Not a valid example number for the *from coefficients* examples.')
end
end