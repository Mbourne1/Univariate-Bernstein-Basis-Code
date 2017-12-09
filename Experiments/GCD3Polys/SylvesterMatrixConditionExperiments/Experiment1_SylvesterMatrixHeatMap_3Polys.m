close all;
clc;


m = 29;
n = 19;
o = 18;
k = 1;


bool_reorder_polynomials = false;
SylvesterMatrixHeatMap_3Polys(m, n, o, k, bool_reorder_polynomials);

m = 19;
n = 29;
o = 18;
k = 1;


bool_reorder_polynomials = false;
SylvesterMatrixHeatMap_3Polys(m, n, o, k, bool_reorder_polynomials);


m = 18;
n = 19;
o = 29;
k = 1;


bool_reorder_polynomials = false;
SylvesterMatrixHeatMap_3Polys(m, n, o, k, bool_reorder_polynomials);

%bool_reorder_polynomials = true;
%SylvesterMatrixHeatMap_3Polys(m, n, o, k, bool_reorder_polynomials);