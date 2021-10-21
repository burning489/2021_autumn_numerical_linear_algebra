function L = func_cholesky(A)
% FUNC_CHOLESKY Cholesky factorization LL^{T}.
% --------------------
% Usage
% L = FUNC_CHOLESKY(A) returns lower triangular matrix L, such that 
% A = L*L'. 
% --------------------
% Input
% A: (n,n) double
%    S.P.D matrix to factorize
% --------------------
% Output
% L: (n,n) double
%    lower triangular matrix

% shape assertions
[m,n] = size(A);
assert(m == n);

A(1,1) = sqrt(A(1,1));
A(2:n,1) = A(2:n,1)/A(1,1);
for j = 2:n
    A(j:n,j) = A(j:n,j)-A(j:n,1:j-1)*A(j,1:j-1)';
    A(j,j) = sqrt(A(j,j));
    A(j+1:n,j) = A(j+1:n,j)/A(j,j);
end
L = tril(A);
end
