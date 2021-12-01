function [x,exitflag,output] = func_conjugate_gradient(A,b,x0,options)
% FUNC_CONJUGATE_GRADIENT Conjugate gradient method to solve SPD problem: Ax = b
% --------------------
% Usage
% [X,EXITFLAT,OUTPUT] = FUNC_CONJUGATE_GRADIENT(A,B,X0) 
% [X,EXITFLAT,OUTPUT] = FUNC_CONJUGATE_GRADIENT(A,B,X0,OPTIONS) 
% returns solution X and OUTPUT information containing number of iterations
% and norms of residuals at each iteration
% --------------------
% Input
% A: (n,n) double
%    SPD matrix
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
%    cg solution
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
    errID = "FUNC_CONJUGATE_GRADIENT:NotEnoughInput";
    msgtext = "func_conjugate_gradient should at least receive A, b, and x0";
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

r_norms = [];
b_norm = norm(b,2);
x = x0;
r = b-A*x;
rho = r'*r;
for iter = 1:max_iter
    % compute direction p
    if (iter==1)
        p = r;
    else
        beta = rho/rho_old;
        p = r+beta*p;
    end
    w = A*p;
    % compute step size
    alpha = rho/(p'*w);
    % one step update
    x = x+alpha*p;
    % compute residual
    r = r-alpha*w;
    rho_old = rho;
    rho = r'*r;
    r_norm = sqrt(rho);
    r_norms = [r_norms,r_norm];
    % check stopping criteria
    if (r_norm<tol*b_norm)
        exitflag = 1;
        output.iter = iter;
        output.r_norms = r_norms;
        return;
    end
end
exitflag = 0;
output.iter = iter;
output.r_norms = r_norms;
end
