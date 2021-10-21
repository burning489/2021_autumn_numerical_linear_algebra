function x = func_ldlt_solver(A,b)
% FUNC_LDLTSOLVER Solve linear equations with LDL^{T} factorization.
% --------------------
% Usage
% x = FUNC_LDLTSOLVER(A,b) return the solution of A*x=b, where A is S.P.D.
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

% LDL^{T} factorization
[L,D] = func_ldlt(A);
% forward substitution to solve L*y=b
y = func_forward(L,b);
% backward substitution to solve D*L'*x=y
y = y./diag(D);
x = func_backward(L',y);
end
