% -------------------------------------------------
%
%
% -------------------------------------------------

clear, clc

n = 200;
A = rand(n,n);
b = rand(n,1);

[L1, U1] = func_lu(A);
err_relative1 = norm(L1*U1-A)/norm(A);

[L2,U2,P] = func_partial_piv_lu(A);
err_relative2 = norm(L2*U2-P*A)/norm(A);

[L3,U3,P,Q] = func_full_piv_lu(A);
err_relative3 = norm(L3*U3-P*A*Q)/norm(A);

x1 = func_lu_solver(A,b,"raw");
x2 = func_lu_solver(A,b,"partial");
x3 = func_lu_solver(A,b,"full");
err_relative4 = norm(A*x1-b)/norm(b);
err_relative5 = norm(A*x2-b)/norm(b);
err_relative6 = norm(A*x3-b)/norm(b);

fprintf("%-60s%f\t%f\n", "relative error of LU factorization without pivoting:", err_relative1, err_relative4);
fprintf("%-60s%f\t%f\n", "relative error of LU factorization with partial pivoting:", err_relative2, err_relative5);
fprintf("%-60s%f\t%f\n", "relative error of raw LU factorization with full pivoting:", err_relative3, err_relative6);