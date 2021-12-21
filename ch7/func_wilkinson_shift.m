function [T] = func_wilkinson_shift(T)
% FUNC_WILKINSON_SHIFT QR method with Wilkinson shift. see p.209
% --------------------
% Usage
% [T] = FUNC_WILKINSON_SHIFT(T)

n = max(size(T));
d = (T(n-1, n-1) - T(n, n))/2;
mu = T(n,n) + d-sign(d) * sqrt(d * d+T(n-1, n) * T(n, n-1));
x = T(1, 1) - mu;
z = T(2, 1);

[c, s] = func_givens(x, z);
G = [c, s; -s, c];
if n == 2
    T(1:2, 1:2) = G' * T(1:2, 1:2);
    T(1:2, 1:2) = T(1:2, 1:2) * G;
else
    T(1:2, 1:3) = G' * T(1:2, 1:3);
    T(1:3, 1:2) = T(1:3, 1:2) * G;
end

for k = 2:n-1
        x = T(k, k-1);
        z = T(k+1, k-1);
        [c, s] = func_givens(x, z);
        G = [c, s; -s, c];
        if k < n-1
            T(k:k+1, k-1:k+2) = G' * T(k:k+1, k-1: k+2);
            T(k-1:k+2, k:k+1) = T(k-1:k+2, k:k+1) * G;
        else
            T(n-1:n,n-2:n) = G' * T(n-1:n, n-2:n);
            T(n-2:n,n-1:n) = T(n-2:n,n-1:n) * G;
        end
end
end