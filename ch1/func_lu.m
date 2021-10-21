function [L,U] = func_lu(A)
% FUNC_LU LU factorization.
% Usage
% --------------------
% [L,U] = FUNC_LU(A) returns the LU of A without pivoting.
% Input
% --------------------
% A: (n,n) double
%    matrix to be LU factorized, should be nonsingual
% Output
% --------------------
% L: (n,n) double
%    unit lower triangular matrix
% U: (n,n) double
%    upper triangular matrix

% shape assertions
[m,n] = size(A);
assert(m == n);

for i = 1:n-1
    if (abs(A(i,i)) < eps)
        % raise error if upper left element is close to 0
        error("LU factorization failed.\nInput matrix is singular!");
    else
        A(i+1:n,i) = A(i+1:n,i)/A(i,i);
        A(i+1:n,i+1:n)=A(i+1:n,i+1:n) - A(i+1:n,i)*A(i,i+1:n);
    end
end

L = tril(A,-1) + eye(n);
U = triu(A);
end
