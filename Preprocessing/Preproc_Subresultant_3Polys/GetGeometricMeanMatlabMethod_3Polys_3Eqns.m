function lambda = GetGeometricMeanMatlabMethod_3PolySubresultant(fx, n_k, o_k)

% Build the partition of the Sylvester matrix
C_f1 = BuildSubresultant_Partition_2Polys(fx, n_k);

C_f2 = BuildSubresultant_Partition_2Polys(fx, o_k);

% Get geometric mean of non-zero entries

SetOfEntries = ...
    [...
        abs(C_f1(C_f1~=0))
        abs(C_f2(C_f2~=0))
    ];

lambda = geomean(SetOfEntries);

end