function [u,lambda,exitflag,output] = func_inverse_power(A,u,options)
% FUNC_INVERSE_POWER Inverse power method to compute largest eigenpair
% --------------------
% Usage
% [U,LAMBDA,EXITFLAT,ITER] = FUNC_INVERSE_POWER(A,U)
% [U,LAMBDA,EXITFLAT,ITER] = FUNC_INVERSE_POWER(A,U,OPTIONS) 
% returns eigenpari, convergence indicator EXITFLAG and OUTPUT information containing number of iterations
% and norms of residuals at each iteration
% --------------------
% Input
% A: (n,n) double
%    matrix
% u: (n,1) double
%    initial guess of eigenvector
% options: struct
%   max_iter: integer
%             maximum number of itertions allowed, default=1e3
%   tol: double
%        tolerance for stopping criteria, default=1e-6   
% --------------------
% Output
% u: (n,1) double
%    approximated eigenvector
% lambda: double
%         approximated eigenvalue
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
if nargin < 2
    errID = "FUNC_INVERSE_POWER:NotEnoughInput";
    msgtext = "func_inverse_power should at least receive A and u";
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
assert((m==n)&&(m==size(u,1))&&(size(u,2)==1));

if isempty(u)
    u = randn(n,1);
end
r_norms = [];
u = u/norm(u);

% estimate of 2-norm of A for convergence test
normA_est = sqrt(norm(A,1) * norm(A,inf));

for iter = 1:max_iter
    addpath(genpath("../../"));
    % solve uk for A*uk = ukm1
    u = func_lu_solver(A,u);
    u = u/norm(u);
    % compute Rayleigh quotient
    Au = A*u;
    lambda = u'*Au;
    % compute residual and check versus tolerance
    r = Au - u*lambda;
    r_norm = norm(r);
    r_norms = [r_norms, r_norm];
    if r_norm < tol*normA_est
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