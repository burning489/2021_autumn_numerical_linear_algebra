function x = func_backward(U,y)
% FUNC_BACKWARD Backward method for upper triangular matrix.
% --------------------
% Usage
%   [X] = FUNC_BACKWARD(U,Y) returns the solution X for U*X=Y.
% --------------------
% Input
% U: (n,n) double
%    upper triangular matrix
% y: (n,1) double
%    right hand side vector
% --------------------
% Output
% x: (n,1) double
%    solution to U*x=y

% shape assertions
[m,n] = size(U);
assert(m == n && n == size(y,1));

x = y;
for j = n:-1:2
    x(j) = x(j)/U(j,j);
    % backward substitution
    x(1:j-1) = x(1:j-1) - x(j)*U(1:j-1,j);
end

x(1) = x(1)/U(1,1);
end
