%this function fills the matrix of the type [0 rot ;rot 0] with
%symetrically computed derivative
function S = FillA_S_block(BC1)
global nx ny nz;

n=nx*ny*nz;
fiz = fi(nz);

% if strcmp(BC1, 'PEC')
%     fiz(1,1) = 0;
% end

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


%FORWARD AND DIAG MATRIX
function Dout = D(n)
    Dout = speye(n);
end
function Dpout = Dp(n)
    Dpout = sparse(1:n-1,2:n,ones(1,n-1),n,n);
    Dpout(n, 1) = 1;
end
function T = nad_diagonalni(n)
    T = (2.0/3.0)*sparse(1:n-1,2:n, ones(1,n-1),n,n) - (1.0/12.0) * sparse(1:n-2,3:n, ones(1,n-2),n,n);
    T(1, n-1) = 1.0/12.0;
    T(2, n) = 1.0/12.0;
    T(1, n) = -2.0/3.0;
end
function fiout = fi(n)
    %fiout = nad_diagonalni(n) - nad_diagonalni(n)' ;
    fiout = 0.5*(Dp(n) - Dp(n)'); %symetric derivative of the first order with periodic boundary conditions
end
%CHANGE OF INDEX
function f = dind(di,dj,dk,dcomp)
    f = 3*di + 3*nx*dj + 3*nx*ny*dk + dcomp;
end
end
