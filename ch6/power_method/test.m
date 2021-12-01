%% simple correctness test
clc, clear

% generate SPD matrix A
n = 5;
rng(42);
A = rand(n,n);
A = A+A';

options.max_iter = 1e3;
options.tol = 1e-10;

u0 = ones(n,1);

%% power method
[x1,lambda1] = func_power(A,u0,options);
[eigvec,eigval] = eig(A);

fprintf([repmat('-', 1, 40), '\n']);
fprintf("power method:\n")
fprintf('max eigen value = %6.4f\n',lambda1);
fprintf('corresponding eigenvector = \n');
fprintf("%6.4f  ", x1);
fprintf("\n");

fprintf([repmat('-', 1, 40), '\n']);
fprintf("MATLAB build-in eig:\n")
fprintf('max eigen value = %6.4f\n', eigval(end));
fprintf('corresponding eigenvectors = \n');
fprintf("%6.4f  ", eigvec(:,end));
fprintf("\n");
fprintf([repmat('-', 1, 40), '\n']);

%% shifted inverse power method
lambda = eigval(3,3);

x = func_shifted_inverse_power(A,u0,lambda,options);
fprintf('shifted inverse power method:\n');

fprintf('A*x = \n');
fprintf("%6.4f  ", A*x);
fprintf("\n");
fprintf('lambda*x = \n');
fprintf("%6.4f  ", A*x);
fprintf("\n");
err = norm(A*x-lambda*x);
fprintf('||A*x-lambda*x|| = %6.4f\n',err);
fprintf([repmat('-', 1, 40), '\n']);