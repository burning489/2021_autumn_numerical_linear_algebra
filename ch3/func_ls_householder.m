function x = func_ls_householder(A,b)
% FUNC_LS_HOUSEHOLDER Solve least square problem with Householder
% transformation.
% --------------------
% Usage
% X = FUNC_LS_HOUSEHOLDER(A,B) returns solution to least square problem
% min |A*x-b|
% --------------------
% Input
% A: (m,n) double
%    matrix
% b: (m,1) double
%    right hand side vector
% --------------------
% Output
% x: (n,1) double
%    solution to least square problem

% shape assertions
[m,n] = size(A);
assert( m >= n && m == length(b));

% qr
[A,d] = func_qr_householder(A);

% Q'b
for j = 1:n
    vb = b(j)+dot(A(j+1:end,j),b(j+1:end));
    beta = d(j);
    b(j) = b(j)-beta*vb;
    b(j+1:end) = b(j+1:end)-beta*vb*A(j+1:end,j);
end

% backward substitution: Rx=Q'b
for j = n:-1:2
    b(j) = b(j)/A(j,j);
    b(1:j-1) = b(1:j-1)-b(j)*A(1:j-1,j);
end
b(1) = b(1)/A(1,1);
x = b(1:n);
end
