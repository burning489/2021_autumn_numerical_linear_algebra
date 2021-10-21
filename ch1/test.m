clear, clc;

%% func_forward
clear
n = 100;
A = rand(n,n)+0.1*n*diag(ones(n,1));
L = tril(A);
b = rand(n,1);
tol = 1.0e-12;

fprintf("test: func_forward\n");
x = func_forward(L,b);
x_star = L\b;
err = norm(x-x_star);
assert(err<tol);
fprintf("\terr = %6.4e\n",err);

%% func_backword
clear 
n = 100;
A = rand(n,n)+0.1*n*diag(ones(n,1));
U = triu(A);
b = rand(n,1);
tol = 1.0e-12;

fprintf("test: func_backword\n");
x = func_backward(U,b);
x_star = U\b;
err = norm(x-x_star);
assert(err<tol);
fprintf("\terr = %6.4e\n",err);

%% func_lu
clear
n = 5;
A = rand(n,n)+0.1*n*diag(ones(n,1));
tol = 1.0e-12;

fprintf("test: func_lu\n");
[L,U] = func_lu(A);
LU = L*U;
err = norm(A(:)-LU(:),inf);
assert(err<tol);
fprintf("\terr = %6.4e\n",err);

%% func_partialpivlu
clear
n = 100;
A = rand(n,n);
tol = 1.0e-12;

fprintf("test: func_partialpivlu\n");
[L,U,P] = func_partial_piv_lu(A);
PA = P*A;
LU = L*U;
err = norm(PA(:)-LU(:),inf);
assert(err<tol);
fprintf("\terr = %6.4e\n",err);

%% func_completepivlu
clear
n = 100;
A = rand(n,n);
tol = 1.0e-12;

fprintf("test: func_completepivlu\n");
[L,U,P,Q] = func_full_piv_lu(A);
PAQ = P*A*Q;
LU = L*U;
err = norm(PAQ(:)-LU(:),inf);
assert(err<tol);
fprintf("\terr = %6.4e\n",err);

%% func_lusolver
clear
n = 100;
A = rand(n,n)+0.1*n*diag(ones(n,1));
b = rand(n,1);
tol = 1.0e-12;

fprintf("test: func_lusolver\n");
x_star = A\b;
x = func_lu_solver(A,b,"raw");
err = norm(x-x_star);
assert(err<tol);
fprintf("\terr = %6.4e\n",err);
x = func_lu_solver(A,b,"partial");
err = norm(x-x_star);
assert(err<tol);
fprintf("\terr = %6.4e\n",err);
x = func_lu_solver(A,b,"full");
err = norm(x-x_star);
assert(err<tol);
fprintf("\terr = %6.4e\n",err);

%% func_cholesky
clear
n = 100;
A = rand(n,n)+0.1*n*diag(ones(n,1));
AS = A+A';
tol = 1.0e-12;

fprintf("test: func_cholesky\n");
L = func_cholesky(AS);
LLT = L*L';
err = norm(LLT(:)-AS(:),inf);
assert(err<tol);
fprintf("\terr = %6.4e\n",err);

%% func_choleskysolver
clear
n = 100;
A = rand(n,n)+0.1*n*diag(ones(n,1));
AS = A+A';
b = rand(n,1);
tol = 1.0e-12;

fprintf("test: func_choleskysolver\n");
x_star = AS\b;
x = func_cholesky_solver(AS,b);
err = norm(x-x_star);
assert(err<tol);
fprintf("\terr = %6.4e\n",err);


%% func_ldlt
clear
n = 100;
A = rand(n,n)+0.1*n*diag(ones(n,1));
AS = A+A';
tol = 1.0e-12;

fprintf("test: func_ldlt\n");
[L,D] = func_ldlt(AS);
LDLT = L*D*L';
err = norm(LDLT(:)-AS(:),inf);
assert(err<tol);
fprintf("\terr = %6.4e\n",err);

%% func_ldltsolver
clear
n = 100;
A = rand(n,n)+0.1*n*diag(ones(n,1));
AS = A+A';
b = rand(n,1);
tol = 1.0e-12;

fprintf("test: func_ldltsolver\n");
x_star = AS\b;
x = func_ldlt_solver(AS,b);
err = norm(x-x_star);
assert(err<tol);
fprintf("\terr = %6.4e\n",err);
