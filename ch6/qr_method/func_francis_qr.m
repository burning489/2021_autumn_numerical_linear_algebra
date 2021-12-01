function H = func_francis_qr(H)
% FUNC_FRANCIS_QR Implicit double shifted QR itreration, see p.193
n = size(H,1);
m = n-1;
s = H(m,m)+H(n,n);
t = H(m,m)*H(n,n)-H(m,n)*H(n,m);
x = H(1,1)*H(1,1)+H(1,2)*H(2,1)-s*H(1,1)+t;
y = H(2,1)*H(1,1)+H(2,2)*H(2,1)-s*H(2,1);
z = H(3,2)*H(2,1);
for k = 0:n-3
    addpath(genpath("../../"));
    [v,b] = func_householder([x,y,z]');
    q = max(1,k);
    H(k+1:k+3,q:n) = H(k+1:k+3,q:n)-b*v*(v'*H(k+1:k+3,q:n));
    r = min(k+4,n);
    H(1:r,k+1:k+3) = H(1:r,k+1:k+3)-b*(H(1:r,k+1:k+3)*v)*v';
    x = H(k+2,k+1);
    y = H(k+3,k+1);
    if k<n-3
        z = H(k+4,k+1);
    end
end
[v,b] = func_householder([x,y]');
H(n-1:n,n-2:n) = H(n-1:n,n-2:n)-b*v*(v'*H(n-1:n,n-2:n));
H(1:n,n-1:n) = H(1:n,n-1:n)-b*(H(1:n,n-1:n)*v)*v';
end
