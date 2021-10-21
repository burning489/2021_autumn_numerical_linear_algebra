function [L,U,P,Q] = func_full_piv_lu(A)
% FUNC_FULL_PIV_LU LU factorization with full pivoting.
% Usage
% --------------------
% [L,U,P,Q] = FUNC_COMPLETEPIVLU(A) returns upper triangular matrix U, 
% lower triangular matrix L and row/column permutation matrix P and Q 
% such that P*A*Q=L*U.
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
% P: (1,n) double
%    row permutation matrix
% Q: (1,n) double
%    column permutation matrix

% shape assertions
[m,n] = size(A);
assert(m == n);

P = eye(n);
Q = eye(n);
for i = 1:n-1
    % find element with largest absolute value in A(i:n,i:n) and its index
    [~,row_idx,col_idx] = find_max_element(abs(A(i:n,i:n)));
    % get index in A
    row_idx = row_idx + (i-1);
    col_idx = col_idx + (i-1);
    if (abs(A(row_idx,col_idx)) < eps)
        % raise error if column is 0
        error("PLU factorization failed.\nInput matrix is singular!");
    else
        % swap row and column and record permutation if necessary
        if (row_idx ~= i || col_idx ~= i)
            A([i,row_idx],:) = A([row_idx,i],:);
            P([i,row_idx],:) = P([row_idx,i],:);
            A(:,[i,col_idx]) = A(:,[col_idx,i]);
            Q(:,[i,col_idx]) = Q(:,[col_idx,i]);
        end
        A(i+1:n,i) = A(i+1:n,i)/A(i,i);
        A(i+1:n,i+1:n)=A(i+1:n,i+1:n)-A(i+1:n,i)*A(i,i+1:n);
    end
end

L = tril(A,-1) + eye(n);
U = triu(A);
end

% get largest element and its index in matrix A 
function [max_val,row_idx,col_idx] = find_max_element(A)
    [col_max,rows] = max(A);
    [max_val,col_idx] = max(col_max);
    row_idx = rows(col_idx);
end
