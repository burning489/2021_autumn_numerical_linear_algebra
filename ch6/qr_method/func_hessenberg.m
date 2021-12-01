function A = func_hessenberg(A)
% FUNC_HESSENBERG Hessenberg decomposition with Householder transform, see p.181
n = size(A,1);
for k = 1:n-2
    addpath(genpath("../../"));
    [v,beta] = func_householder(A(k+1:n,k));
    A(k+1:n,k:n) = A(k+1:n,k:n)-beta*v*(v'*A(k+1:n,k:n));
    A(1:n,k+1:n) = A(1:n,k+1:n)-beta*(A(1:n,k+1:n)*v)*v';
end
end
