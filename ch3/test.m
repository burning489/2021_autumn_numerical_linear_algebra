clear, clc, close all;

%% Householder
clear
n = 10;
x = rand(n,1);
[v,beta] = func_householder(x);
v(1) = 1.0;
H = eye(n)-beta*(v*v');
Hx = H*x;
err = norm(Hx(2:n));
fprintf("test: func_householder\n");
fprintf("\terr = %6.4e\n",err);
tol = 1.0e-12;
assert(err < tol)

%% QR-GramSchmidt
clear
m = 10;
n = 4;
A = rand(m,n);
[Q,R] = func_qr_gramschmidt(A);
QR = Q*R;
err = norm(QR(:)-A(:),inf);
fprintf("test: func_qr_gramschmidt\n");
fprintf("\terr = %6.4e\n",err);
tol = 1.0e-12;
assert(err < tol);

%% QR-Householder
clear
m = 10;
n = 4;
A = rand(m,n);
[AA,d] = func_qr_householder(A);

v1 = [1.0;AA(2:end,1)];
H1 = eye(m)-d(1)*(v1*v1');
A1 = H1*A;
v2 = [0.0;1.0;AA(3:end,2)];
H2 = eye(m)-d(2)*(v2*v2');
A2 = H2*A1;
v3 = [0.0;0.0;1.0;AA(4:end,3)];
H3 = eye(m)-d(3)*(v3*v3');
A3 = H3*A2;
v4 = [0.0;0.0;0.0;1.0;AA(5:end,4)];
H4 = eye(m)-d(4)*(v4*v4');
A4 = H4*A3;

% A=Q*R
Q = H1*H2*H3*H4;
R = triu(AA);
QR = Q*R;
err = norm(QR(:)-A(:),inf);
fprintf("test: func_qr_householder\n");
fprintf("\terr = %6.4e\n",err);
tol = 1.0e-12;
assert(err < tol);

%% least square
clear
m = 100;
n = 4;
A = rand(m,n);
b = rand(m,1);
x_house = func_ls_householder(A,b);
x_schmidt = func_ls_gramschmidt(A,b);
xx = A\b;

tol = 1.0e-12;
err_house = norm(x_house-xx,inf);
fprintf("test: func_ls_householder\n");
fprintf("\terr = %6.4e\n",err_house);
assert(err_house < tol);

err_schmidt = norm(x_schmidt-xx,inf);
fprintf("test: func_ls_gramschmidt\n");
fprintf("\terr = %6.4e\n",err_schmidt);
assert(err_schmidt < tol);

