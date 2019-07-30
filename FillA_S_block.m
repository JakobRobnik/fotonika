%this function fills the matrix of the type [0 rot ;rot 0] with
%symetrically computed derivative
function S = FillA_S_block(BC1)
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


rotor = [sparse(n,n) -Fz Fy; Fz sparse(n,n) -Fx; -Fy Fx sparse(n,n)];
S = [sparse(3*n, 3*n) rotor; rotor sparse(3*n, 3*n)];

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
    Dpout(n, 1) = 1;
end
function fiout = fi(n)
    fiout = 0.5*(Dp(n) - Dp(n)'); %symetric derivative with periodic boundary conditions
end
%CHANGE OF INDEX
function f = dind(di,dj,dk,dcomp)
    f = 3*di + 3*nx*dj + 3*nx*ny*dk + dcomp;
end
end
