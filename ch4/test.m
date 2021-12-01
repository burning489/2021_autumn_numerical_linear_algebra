%% simple correctness test
clc, clear, close all
rng(42);
n = 10;
A = rand(n,n);
A = A+A'+10.0*eye(n);
b = rand(n,1);
x0 = zeros(n,1);

max_iter = 1e3;
tol = 1e-10;

options.max_iter = max_iter;
options.tol = tol;

% Jacobi 
[x_jac, exitflag_jac, output_jac] = func_jacobi(A,b,x0,options);

% Gauss-Seidel
[x_gs, exitflag_gs, output_gs] = func_gauss_seidel(A,b,x0,options);


% Display results
fprintf([repmat('-', 1, 40), '\n']);
fprintf("Jacobi:\n");
fprintf("Number of iterations: %d\nError in 2-norm: %.12f\n", output_jac.iter, norm(output_jac.r_norms(end)));
fprintf([repmat('-', 1, 40), '\n']);
fprintf("Gauss Seidel:\n");
fprintf("Number of iterations: %d\nError in 2-norm: %.12f\n", output_gs.iter, norm(output_gs.r_norms(end)));
fprintf([repmat('-', 1, 40), '\n']);

nexttile;
plot(1:output_jac.iter, output_jac.r_norms);
title("Jacobi iterative method")
xlabel("iter")
ylabel("2-norm of residuals")
nexttile;
plot(1:output_gs.iter, output_gs.r_norms);
title("Gauss-Seidel iterative method")
xlabel("iter")
ylabel("2-norm of residuals")
