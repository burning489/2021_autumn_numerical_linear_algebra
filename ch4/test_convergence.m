%% a simple numerical verification for exercise 1 in chapter 4, p.134
clc, clear, close all
A1 = [2 -1 1 ; 1 1 1; 1 1 -2];
A2 = [1 2 -2; 1 1 1; 2 2 1];
b = [-1 2 3]';
x0 = [0 0 0]';
n = 3;

max_iter = 1e3;
tol = 1e-10;

options.max_iter = max_iter;
options.tol = tol;

[x_jac, exitflag_jac, output_jac] = func_jacobi(A1,b,x0,options);
[x_gs, exitflag_gs, output_gs] = func_gauss_seidel(A1,b,x0,options);

fprintf([repmat('_', 1, 40),'\n']);
fprintf("For test matrix A1:\n");
disp(A1);
if exitflag_jac
    fprintf("Jacobi method converged\n");
else
    fprintf("Jacobi method did not converge\n");
end
if exitflag_gs
    fprintf("Gauss-Seidel method converged\n");
else
    fprintf("Gauss-Seidel method did not converge\n");
end


[x_jac, exitflag_jac, output_jac] = func_jacobi(A2,b,x0,options);
[x_gs, exitflag_gs, output_gs] = func_gauss_seidel(A2,b,x0,options);

fprintf([repmat('_', 1, 40),'\n']);
fprintf("For test matrix A2:\n");
disp(A2);
if exitflag_jac
    fprintf("Jacobi method converged\n");
else
    fprintf("Jacobi method did not converge\n");
end
if exitflag_gs
    fprintf("Gauss-Seidel method converged\n");
else
    fprintf("Gauss-Seidel method did not converge\n");
end

fprintf([repmat('_', 1, 40),'\n']);
