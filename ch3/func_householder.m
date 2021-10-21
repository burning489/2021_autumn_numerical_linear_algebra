function [x,beta] = func_householder(x)
% FUNC_HOUSEHOLDER Householder transformation on vector x.
% --------------------
% Usage
% [X,BETA] = FUNC_HOUSEHOLDER(X) returns x and beta corresponding to
% some householder transformation H = I - beta*x*x'.
% --------------------
% Input
% x: (n,1) double
%    vector to be transformed
% --------------------
% Output
% x: (n,1) double
% beta: double
%       such that H = I - beta*x*x'

n = length(x);
eta = max(x);
x = x/eta;
sigma = dot(x(2:end),x(2:end));
if sigma < eps
    beta = 0.0;
else
    alpha = sqrt(x(1)^2 + sigma);
    if x(1) <= 0.0
        x1 = x(1) - alpha;
    else
        x1 = -sigma/(x(1)+alpha);
    end
    beta = 2*x1^2/(sigma+x1^2);
    x(2:n) = x(2:n)/x1;
    x(1) = 1;
end
