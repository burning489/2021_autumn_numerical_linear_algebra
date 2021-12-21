function [D] = func_sym_qr(A, tol)
% FUNC_SYM_QR QR methods to compute eigenvalues of symmetrix matrix A.
% --------------------
% Usage
% [D] = FUNC_SYM_JACOBI(A)
% --------------------
% Input
% A: (n,n) double
%    symmetirx matrix
% tol: double
%      tolerance
% --------------------
% Output
% D: (n,n) double
%    eigenvalues

n = max(size(A));
A = func_householder_tridiag(A);
A = tril(A);
A = A + diag(diag(A, -1), 1);
D = A;
q = 0;

while q < n
    for i = 1:n - 1
        if abs(D(i+1, i)) <= tol * (abs(D(i, i))+abs(D(i+1, i+1)))
            D(i+1, i) = 0;
        end
        [p, q] = select_index(D);
        if q < n
            D(p+1:n - q, p+1:n - q) = func_wilkinson_shift(D(p+1:n - q, p+1:n - q));
        end
    end
end
end

function [p,q] = select_index(D)
n = max(size(D));
q = 0;
n_p_q = 0;
flag = 0;

for i=n-1:-1:1
    if flag == 0 && D(i+1,i) ~= 0       
        flag = 1;
        n_p_q = 1;
    end
    
    if flag == 1 && D(i+1,i) == 0
        break;
    end
    if flag == 0
        q = q+1;
    end
    if flag==1
        n_p_q = n_p_q+1;
    end
end
if q == n-1 && n_p_q == 0
    q = n;
end
p = n-n_p_q - q;
end