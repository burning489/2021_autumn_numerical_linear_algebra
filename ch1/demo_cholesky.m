% Demonstrate correctness of Cholesky factorization and LDL^{T} factorization
% on random 200*200 S.P.D matrix A and solve equation A*x=b. Compute the relative error.
clear, clc

% initialize matrix A and lhs vector b
n = 200;
A = rand(n,n);
A = A + A' + diag(0.1*n*ones(1,n));
b = rand(n,1);

L1 = func_cholesky(A);
err_relative1 = norm(L1*L1'-A)/norm(A);

[L2, D2] = func_ldlt(A);
err_relative2 = norm(L2*D2*L2'-A)/norm(A);

x1 = func_cholesky_solver(A,b);
x2 = func_ldlt_solver(A,b);
err_relative3 = norm(A*x1-b)/norm(b);
err_relative4 = norm(A*x2-b)/norm(b);

fprintf("%-60s%f\t%f\n", "relative error of Cholesky factorization:", err_relative1, err_relative3);
fprintf("%-60s%f\t%f\n", "relative error of LDL^{T} factorization:", err_relative2, err_relative4);
