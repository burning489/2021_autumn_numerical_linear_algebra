function [x,exitflag,output] = func_gauss_seidel(A,b,x0,options)
% FUNC_GAUSS_SEIDEL Gauss-Seidel iterative method to solve Ax = b
% --------------------
% Usage
% [X,EXITFLAT,OUTPUT] = FUNC_GAUSS_SEIDEL(A,B,X0) 
% [X,EXITFLAT,OUTPUT] = FUNC_GAUSS_SEIDEL(A,B,X0,OPTIONS) 
% returns solution X and OUTPUT information containing number of iterations
% and norms of residuals at each iteration
% --------------------
% Input
% A: (n,n) double
%    matrix to solve
% b: (n,1) double
%    right hand side vector
% x0: (n,1) double
%     initial guess of solution
% options: struct
%   max_iter: integer
%             maximum number of itertions allowed, default=1e3
%   tol: double
%        tolerance for stopping criteria, default=1e-6   
% --------------------
% Output
% x: (n,1) double
%    gauss-seidel iterative solution
% exitflag: integer
%           0 maximum number of iterations reached
%           1 perfectly worked 
% output: struct
%   iter: integer
%         number of iterations
%   r_norms: (iter,1) double
%            2-norm of residuals at each iteration

% prepare parameters
if ~exist('options','var')
    options = [];
end
if nargin < 3
    errID = "FUNC_GAUSS_SEIDEL:NotEnoughInput";
    msgtext = "func_gauss_seidel should at least receive A, b, and x0";
    ME = MException(errID,msgtext);
    throw(ME);
end
if ~isfield(options,'max_iter')
    max_iter = 1e3;
else
    max_iter = options.max_iter;
end
if ~isfield(options,'tol')
    tol = 1e-6;
else
    tol = options.tol;
end

% shape assertions
[m,n] = size(A);
assert((m==n)&&(m==size(b,1))&&(size(b,2)==1));

% Get estimate of 2-norm of A for convergence test
normA_est = sqrt(norm(A,1) * norm(A,inf));

r_norms = [];
xk = x0;
xkp1 = zeros(size(xk));
for iter = 1:max_iter
    for i=1:n
        sigma = 0;
        for j = 1:i-1
            sigma = sigma + A(i,j)*xkp1(j);
        end
        for j = i+1:n
            sigma = sigma + A(i,j)*xk(j);
        end
        xkp1(i) = (b(i)-sigma)/A(i,i);
    end
    delta = norm(xkp1-xk);
    xk = xkp1;
    r_norms = [r_norms,norm(A*xk-b)];
    % check stopping criteria
    if (delta < tol*normA_est)
        x = xk;
        exitflag = 1;
        output.iter = iter;
        output.r_norms = r_norms;
        return
    end
end
x = xk;
exitflag = 0;
output.iter = iter;
output.r_norms = r_norms;
end
