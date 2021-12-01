function H = func_schur(H, u)
% FUNC_SCHUR Implicit QR with Schur decomposition, see p.194
% note u is user-defined epsilon
if ~exist('u','var')
    u = 1e-6;
end
n = size(H,1); 
H = func_hessenberg(H);
m = 0;
while 1
    for k = 1:n-1
        if abs(H(k+1,k))<u*(abs(H(k,k))+abs(H(k+1,k+1)))
            H(k+1,k) = 0.0;
        end
    end
    for k = n-m:-1:2
        if H(k,k-1)==0
            m = m+1;
        elseif k == 2
            m = m+2;
        elseif H(k-1,k-2)==0
            m = m+1;
        else
            break
        end
    end
    l = 0;
    for k = n-m:-1:2
        if abs(H(k,k-1))<u
            l = k-1;
            break
        end
    end
    if m >= n-2
        break
    else
        H(l+1:n-m,l+1:n-m) = func_francis_qr(H(l+1:n-m,l+1:n-m));
    end
end

end