function x = func_forward(L,y)
% FUNC_FORWARD Forward method for lower triangular matrix.
% --------------------
% Usage
% [X] = FUNC_FORWARD(L,Y) returns the solution X for L*X=Y.
% --------------------
% Input
% L: (n,n) double
%    lower triangular matrix
% y: (n,1) double
%    right hand side vector
% --------------------
% Output
% x: (n,1) double
%    solution to L*x=y

% shape assertions
[m,n] = size(L);
assert(m == n && n == size(y,1));

x = y;
for j = 1:n-1
    x(j) = x(j)/L(j,j);
    % forward substitution
    x(j+1:n) = x(j+1:n) - x(j)*L(j+1:n,j);
end

x(n) = x(n)/L(n,n);
end