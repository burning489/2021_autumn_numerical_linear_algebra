clear, clc

n = 10;
rng(42);
A = rand(n,n);
A = A+A';

eigval = func_eigval(A);
[~, eigval_exact] = eig(A);
fprintf("eigenvalues computed with Schur:\n");
fprintf("%6.4f  ", sort(eigval,'descend'));
fprintf("\n")
fprintf("eigenvalues computed with MATLAB build-in eig:\n");
fprintf("%6.4f  ", sort(diag(eigval_exact),'descend'));
fprintf("\n")

