function [Q,R] = func_qr_gramschmidt(A)
% FUNC_QR_GRAMSCHMIDT Gram-Schmidt orthogonalization.
% --------------------
% Usage
% [Q,R] = FUNC_QR_GRAMSCHMIDT(A) returns orthogonal matrix Q such that
% A=Q*R 
% --------------------
% Input
% A: (m,n) double
%    matrix to be orthogonalized, should be overdetermined
% --------------------
% Output
% Q: (m,n) double
%    column-wise orthogonal matrix
% R: (n,n) double
%    lower triangular matrix

% shape assertions
[m,n] = size(A);
assert(m >= n);

Q = A;
R = zeros(n,n);
for k = 1:n
    R(k,k) = norm(Q(:,k));
    Q(:,k) = Q(:,k)/R(k,k);
    for j = k+1:n
        R(k,j) = dot(Q(:,k),Q(:,j));
        Q(:,j) = Q(:,j)-R(k,j)*Q(:,k);
    end
end
end
