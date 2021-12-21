function [D,U] = func_sym_jacobi(A,tol)
% FUNC_SYM_JACOBI Jacobi methods to compute eigenvalues of symmetrix matrix A.
% --------------------
% Usage
% [D] = FUNC_SYM_JACOBI(A) to compute eigenvalues
% [D,U] = FUNC_SYM_JACOBI(A) to compute eigenvalues and eigenvectors
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
% U: (n,n) double
%    corresponding eigenvector

n = max(size(A));
U = eye(n);
l = 0;
br = 1;
while br == 1
    for i = 1:n-1
        for j = i+1:n
            if(abs(A(i,j)) >= tol)
                zeta = (A(j,j)-A(i,i))/(2*A(i,j));
                tau = sign(zeta)/(abs(zeta)+sqrt(1+zeta*zeta));
                c = 1/sqrt(1+tau*tau);
                s = tau*c;
                Urot = eye(n);
                Urot(i,i) = c;
                Urot(j,j) = c;
                Urot(i,j) = s;
                Urot(j,i) = -s;
                A = Urot'*A*Urot;
                U = U*Urot;
                l = 0;
            else
                l = l+1;
                if l == n*(n-1)/2
                    D = A;
                    return;
                end
            end
        end
    end
end
D = A;
end