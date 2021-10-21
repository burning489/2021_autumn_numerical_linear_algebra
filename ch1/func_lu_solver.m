function x = func_lu_solver(A,b,method)
% FUNC_LU_SOLVER Solve linear equations with LU factorization.
% --------------------
% Usage
% x = FUNC_LUSOLVER(A,b,option) return the solution of A*x=b.
% --------------------
% Input
% A: (n,n) double
%    matrix to solve
% b: (n,1) double
%    left hand side vector
% method: string, default="partial", optional
%         method to LU factorization
%         "raw" for LU factorization without pivoting
%         "partial" LU factorization with partial pivoting
%         "full" LU factorization with full pivoting
% --------------------
% Output
% x: (n,1) double
%    solution to A*x=b

% default solve option
if (nargin == 2)
    method = "partial";
end

% shape assertions
[m,n] = size(A);
assert(n == m && n == size(b,1));

switch method
    case "raw"
        % LU factorization without pivoting.
        [L,U] = func_lu(A);
        % forward substitution
        y = func_forward(L,b);
        % backward substitution
        x = func_backward(U,y);
    case "partial"
        % LU factorization with partial pivoting.
        [L,U,P] = func_partial_piv_lu(A);
        % permutate right-hand-side vector with P and forward substitution
        y = func_forward(L,P*b);
        % backward substitution
        x = func_backward(U,y);
    case "full"
        % LU factorization with full pivoting.
        [L,U,P,Q] = func_full_piv_lu(A);
        % permutate right-hand-side vector with P and forward substitution
        y = func_forward(L,P*b);
        % backward substitution
        x = func_backward(U,y);
        % permutate solution x with Q
        x = Q*x;
    otherwise
        error("Invalid sovler option in func_lu_solver!")
end
end
