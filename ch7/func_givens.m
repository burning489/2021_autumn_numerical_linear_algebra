function [c, s] = func_givens(a, b)
% FUNC_GIVENS Compute Givens transformation matrix. see p.90
% --------------------
% Usage
% [C,S] = FUNC_GIVENS(A,B) returns cos and sin in Givens transformation matrix
% --------------------
% Output
% c: double
%    cos
% s: double
%    sin

if b == 0
    c = 1;
    s = 0;
    return;
end
if abs(b) > abs(a)
     tau = -a/b;
     s = 1/sqrt(1 + tau^2);
     c = s * tau;
else 
        tau = -b/a;
        c = 1 / sqrt(1 + tau^2);
        s = c * tau;
end
end