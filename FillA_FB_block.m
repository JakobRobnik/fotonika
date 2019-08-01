%this function fills the matrix
function double_rotor = FillA_FB_block(BC1)
global nx ny nz;

n=nx*ny*nz;
fiz = fi(nz);

if strcmp(BC1, 'PEC')
    fiz(1,1) = 0;
end

Fx = kron(D(nz),kron(D(ny),fi(nx)));
Fy = kron(D(nz),kron(fi(ny),D(nx)));
Fz = kron(fiz,kron(D(ny),D(nx)));

%mirror BC for electric field at the start of domain
% if strcmp(BC1, 'PEC')
% %     Fx(1,1) = 0;
% %     Fy(1,1) = 0;
%     Fz(1,1) = 0;
% end

Bx = -Fx';
By = -Fy';
Bz = -Fz';

F = [sparse(n,n) -Fz Fy; Fz sparse(n,n) -Fx; -Fy Fx sparse(n,n)];
B = [sparse(n,n) -Bz By; Bz sparse(n,n) -Bx; -By Bx sparse(n,n)];

double_rotor = [sparse(3*n, 3*n) F; F sparse(3*n, 3*n)];

%INDEX COUNTER
function f = ind(i,j,k)
    f = (k-1)*nx*ny + (j-1)*nx + (i-1) + 1;
end
%FORWARD AND DIAG MATRIX
function Dout = D(n)
    Dout = speye(n);
end
function Dpout = Dp(n)
    Dpout = sparse(1:n-1,2:n,ones(1,n-1),n,n);
end
function fiout = fi(n)
    fiout = Dp(n) - D(n);
    fiout(n, 1) = 1;
end
%CHANGE OF INDEX
function f = dind(di,dj,dk,dcomp)
    f = 3*di + 3*nx*dj + 3*nx*ny*dk + dcomp;
end
end

