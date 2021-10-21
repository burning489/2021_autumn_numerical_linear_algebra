function [L,U,P] = func_partial_piv_lu(A)
% FUNC_PARTIAL_PIV_LU LU factorization with partial pivoting.
% Usage
% --------------------
% [L,U,P] = FUNC_PLU(A) returns upper triangular matrix U, 
% lower triangular matrix L, and row permutation matrix P
% such that P*A=L*U.
% Input
% --------------------
% A: (n,n) double
%    matrix to be PLU factorized, should be nonsingual
% Output
% --------------------
% L: (n,n) double
%    unit lower triangular matrix
% U: (n,n) double
%    upper triangular matrix
% P: (n,n) double
%    row permutation matrix

% shape assertions
[m,n] = size(A);
assert(m == n);

P = eye(n);
for i = 1:n-1
    % find element with largest absolute value in A(i:n,i) and its row index
    [~,row_idx] = max(abs(A(i:n,i)));
    % get row index in A
    row_idx = row_idx + (i-1);
    if (abs(A(row_idx,i)) < eps)
        % raise error if column is 0
        error("PLU factorization failed.\nInput matrix is singular!");
    else
         % swap row and record permutation if necessary
         if (row_idx ~= i)
            A([i,row_idx],:) = A([row_idx,i],:);
            P([i,row_idx],:) = P([row_idx,i],:);
        end
        A(i+1:n,i) = A(i+1:n,i)/A(i,i);
        A(i+1:n,i+1:n)=A(i+1:n,i+1:n)-A(i+1:n,i)*A(i,i+1:n);
    end
end

L = tril(A,-1) + eye(n);
U = triu(A);
end
