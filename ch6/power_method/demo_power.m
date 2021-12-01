%% compare power method and MATLAB build-in function eig
clear, clc
addpath(genpath("../../"));

rng(42);
n = 6;
A = rand(n,n);
% real symmetric matrix with all real eigenvalues
A = A+A';

options.max_iter = 1e3;
options.tol = 1e-10;

%% MATLAB eig
[~,eig_matlab] = eig(A);

%% power method
A1 = A;
[x1,lambda1] = func_power(A1,ones(n,1),options);
[v1,beta1] = func_householder(x1);
A2 = func_householder_trans(A1,v1,beta1);

A2 = A2(2:end,2:end);
[x2,lambda2] = func_power(A2,ones(n-1,1),options);
[v2,beta2] = func_householder(x2);
A3 = func_householder_trans(A2,v2,beta2);

A3 = A3(2:end,2:end);
[x3,lambda3] = func_power(A3,ones(n-2,1),options);
[v3,beta3] = func_householder(x3);
A4 = func_householder_trans(A3,v3,beta3);

A4 = A4(2:end,2:end);
[x4,lambda4] = func_power(A4,ones(n-3,1),options);
[v4,beta4] = func_householder(x4);
A5 = func_householder_trans(A4,v4,beta4);

A5 = A5(2:end,2:end);
[x5,lambda5] = func_power(A5,ones(n-4,1),options);
[v5,beta5] = func_householder(x5);
A6 = func_householder_trans(A5,v5,beta5);

A6 = A6(2:end,2:end);
[x6,lambda6] = func_power(A6,ones(n-5,1),options);

fprintf("eigenvalues using MATLAB eig:\n")
fprintf("%6.4f  ", sort(diag(eig_matlab), 'descend'));
fprintf("\n")
fprintf("eigenvalues using power method:\n")
eigval = [lambda1, lambda2, lambda3, lambda4, lambda5, lambda6];
fprintf(repmat('%6.4f  ',1,6), sort(eigval, 'descend'));
fprintf("\n")


%% auxiliary functions
function A = func_householder_trans(A,v,beta)
    Av = A*v;
    A = A-beta*v*(v'*A)-beta*Av*v'+beta*beta*(v'*Av)*(v*v');
end
