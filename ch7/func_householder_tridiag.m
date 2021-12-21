function [A] = func_householder_tridiag(A)
% FUNC_HOUSEHOLDER_TRIDIAG Tridiagonalization for symmetric matrix A 
% using Householder Transformation. see p.207
% --------------------
% Usage
% [A] = FUNC_HOUSEHOLDER_TRIDIAG(A) 
% --------------------
% Input
% A: (n,n) double
%    symmetirx matrix
% --------------------
% Output
% A: (n,n) double
%    tridiagonalized A
addpath(genpath("../"));
n = max(size(A));
for k = 1:n - 2
    [v, beta] = func_householder(A(k+1:n, k)); 
    H = eye(n-k) - beta * (v*v');
    A(k+1:n, k) = H * A(k+1:n, k);
    u = beta * A(k+1:n, k+1:n) * v;
    w = u - (beta *u' * v / 2) * v;
    A(k+1, k) = norm(A(k+1:n, k));
    A(k, k+1) = A(k+1, k);
    A(k+1:n, k+1:n) = A(k+1:n, k+1:n) - v * w' - w * v';
end
end
