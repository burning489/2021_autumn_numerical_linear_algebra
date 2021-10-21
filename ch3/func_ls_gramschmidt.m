function x = func_ls_gramschmidt(A,b)
% FUNC_LS_GRAMSCHMIDT Solve least square problem with Gram-Schmidt process.
% --------------------
% Usage
% X = FUNC_LS_GRAMSCHMIDT(A,B) returns solution to least square problem
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
assert(m >= n && m == length(b));

% Gram-Schmidt process
[Q,R] = func_qr_gramschmidt(A);

% backward substitution: Rx=Q'b
b = Q'*b;
for j = n:-1:2
    b(j) = b(j)/R(j,j);
    b(1:j-1) = b(1:j-1)-b(j)*R(1:j-1,j);
end
b(1) = b(1)/R(1,1);
x = b;
end
