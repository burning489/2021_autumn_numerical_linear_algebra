clc, clear
tol = 1e-10;
x = [1, 2, 3, 4];
D = diag(x);
Q = rand(4);
[Q, ~] = qr(Q);
A = Q * D * Q';
eig_qr = func_sym_qr(A, tol);
display(eig_qr);
[eig_jacobi, ~] = func_sym_jacobi(A,tol);
display(eig_jacobi);