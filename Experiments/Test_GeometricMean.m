
fx = [ ...
    1.0654 ;
    57.00 ;
    100003545.1245;
    -1.1546547897654;
    -552.2687654;
    10.15654
    -5753.21546
    121.24654
    654630.35464
    ];

m = GetDegree(fx);


n_k = 20;

global SETTINGS
SETTINGS.SYLVESTER_BUILD_METHOD = 'DTQ';

mean_method = 'Arithmetic Mean';

switch mean_method
    
    case 'Geometric Mean'
        
        lambda = GetGeometricMean(fx,n_k);
        
    case 'Arithmetic Mean'
        
        lambda = GetArithmeticMean(fx,n_k);
        
end

fprintf('Mean : %2.4f \n', lambda)

fx_n = fx ./ lambda;

C1 = BuildSubresultant_Partition_2Polys(fx, n_k);

C1_normalised = BuildSubresultant_Partition_2Polys(fx_n, n_k);

nEntries = (n_k + 1) * (m + n_k + 1);

Entries_C1              = reshape(C1, nEntries,1);
Entries_C1_normalised   = reshape(C1_normalised, nEntries,1);

Entries_C1              = Entries_C1(Entries_C1~=0);
Entries_C1_normalised   = Entries_C1_normalised(Entries_C1_normalised~=0);

figure()
hold on
plot(1:1:length(Entries_C1), log10(abs(Entries_C1)),'-b')
plot(1:1:length(Entries_C1_normalised), log10(abs(Entries_C1_normalised)),'-r')
hold off


delta = max(log10(abs(Entries_C1))) - min(log10(abs(Entries_C1)));
delta_normalised = max(log10(abs(Entries_C1_normalised))) - min(log10(abs(Entries_C1_normalised)));

display(delta)
display(delta_normalised)
display('end')