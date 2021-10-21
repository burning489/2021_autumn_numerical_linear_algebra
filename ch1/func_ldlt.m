function [L,D] = func_ldlt(A)
% FUNC_LDLT LDL^{T} factorization LDL^{T}.
% --------------------
% Usage
% [L,D] = FUNC_LDLT(A) returns unit upper triangular matrix L, 
%   and diagonal matrix D such that A = L*D*L'. 
% --------------------
% Input
% A: (n,n) double
%    S.P.D matrix to factorize
% --------------------
% Output
% L: (n,n) double
%    unit lower triangular matrix
% D: (n,n) double
%    diagonal matrix

% shape assertions
[m,n] = size(A);
assert(m == n);

v = zeros(n,1);
for j = 1:n
    for i = 1:j-1
        v(i) = A(j,i)*A(i,i);
    end
    A(j,j) = A(j,j)-A(j,1:j-1)*v(1:j-1);
    A(j+1:n,j) = A(j+1:n,j)-A(j+1:n,1:j-1)*v(1:j-1);
    A(j+1:n,j) = A(j+1:n,j)/A(j,j);
end

L = tril(A,-1) + eye(n);
D = diag(diag(A));
end
