function [Q,d] = func_qr_householder(A)
% FUNC_QR_HOUSEHOLDER QR factorization using householder transformation.
% --------------------
% Usage
% [A,d] = FUNC_LDLT(A) returns the column-wise orthogonal matrix Q, and d store beta calculated.
% --------------------
% Input
% A: (m,n) double
%    matrix to factorize, should be overdetermined
% --------------------
% Output
% Q: (n,n) double
%    column-wise orthogonal matrix
% d: (n,1) double
%    store beta

% shape assertions
[m,n] = size(A);
assert(m >= n);

Q = A;
d = zeros(n,1);
for j = 1:n
    [v,beta] = func_householder(Q(j:m,j));
    Q(j:m,j:n) = Q(j:m,j:n) - beta*v*v'*Q(j:m,j:n);
    d(j) = beta;
    Q(j+1:m,j) = v(2:end);
end
end
