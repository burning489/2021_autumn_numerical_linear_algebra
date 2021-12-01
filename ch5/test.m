%% simple correctness test
clear, clc
% Generate SPD problem
rng(42);
n = 500;
L = tril(rand(n,n))+n*eye(n);
A = L*L';
b = rand(n,1);
x_exact = A\b;

%% Gradient Descent
x_init = zeros(n,1);
options.max_iter = 100;
options.tol = 1e-10;
[x_gd, flag_gd, output_gd] = func_gradient_decent(A,b,x_init,options); 
err_gd = norm(x_gd-x_exact,2);

%% Conjugate Gradient
x_init = zeros(n,1);
[x_cg, flag_cg, output_cg] = func_conjugate_gradient(A,b,x_init,options); 
err_cg = norm(x_gd-x_exact,2);

% Display results
fprintf([repmat('-', 1, 40), '\n']);
fprintf("Gradient Descent:\n");
fprintf("Number of iterations: %d\nError in 2-norm: %.12f\n", output_gd.iter, err_gd);
fprintf([repmat('-', 1, 40), '\n']);
fprintf("Conjugate Gradient:\n");
fprintf("Number of iterations: %d\nError in 2-norm: %.12f\n", output_cg.iter, err_cg);
fprintf([repmat('-', 1, 40), '\n']);