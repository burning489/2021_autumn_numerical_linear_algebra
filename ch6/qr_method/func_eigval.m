function [eigval] = func_eigval(A)
% FUNC_EIGVAL
% Call FUNC_SCHUR
A = func_schur(A);
eigval = func_schur2eigval(A);
end

function eigval = func_schur2eigval(H)
n = size(H,1); 
eigval = zeros(n,1);
i=1;
while i<=n-1
    if H(i+1,i)==0
        eigval(i) = H(i,i);
        i = i+1;
    else
        b = -H(i,i)-H(i+1,i+1);
        c = H(i,i)*H(i+1,i+1)-H(i,i+1)*H(i+1,i);
        delta = b*b-4*c;
        eigval(i) = 0.5*(-b+sqrt(delta));
        eigval(i+1) = 0.5*(-b-sqrt(delta));
        i = i+2;
    end
end 
if i==n
    eigval(n) = H(n,n);
end
end
