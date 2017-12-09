function lambda = GetGeometricMeanMatlabMethod(fx,n_k)

% Build the partition of the Sylvester matrix
C_f_unproc = BuildSubresultant_Partition_2Polys(fx, n_k);

% Get geometric mean of non-zero entries
lambda = geomean(abs(C_f_unproc(C_f_unproc~=0)));

end