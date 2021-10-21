function x = func_cholesky_solver(A,b)
% FUNC_CHOLESKYSOLVER Solve linear equations with Cholesky factorization.
% --------------------
% Usage
% x = FUNC_CHOLESKYSOLVER(A,b) return the solution of A*x=b, where A is S.P.D.
% --------------------
% Input
% A: (n,n) double
%    S.P.D matrix to solve
% b: (n,1) double
%    right hand side vector
% --------------------
% Output
% x: (n,1) double
%    solution to A*x=y

% shape assertions 
[m,n] = size(A);
assert(n == m && n == size(b,1));

% Cholesky factorization
L = func_cholesky(A);
% forward substitution to solve L*y=b
y = func_forward(L,b);
% backward substitution to solve L'*x=y
x = func_backward(L',y);
end
