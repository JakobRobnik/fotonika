%this function fills the matrix
function double_rotor = FillA_FB_block()
global nx ny nz modetype kvec;

n=nx*ny*nz;

Fx = kron(D(nz),kron(D(ny),fiX(nx,kvec)));
Fy = kron(D(nz),kron(fiY(ny,kvec),D(nx)));
Fz = kron(fiZ(nz,kvec),kron(D(ny),D(nx)));

Bx = -Fx';
By = -Fy';
Bz = -Fz';

if strcmp(modetype, 'ALL') 
    F = [sparse(n,n) -Fz Fy; Fz sparse(n,n) -Fx; -Fy Fx sparse(n,n)];
    B = [sparse(n,n) -Bz By; Bz sparse(n,n) -Bx; -By Bx sparse(n,n)];
elseif strcmp(modetype, 'TM')
    F = [sparse(n,n) -Fz Fy; Fz sparse(n,n) -Fx; sparse(n,n) sparse(n,n) sparse(n,n)];
    B = [sparse(n,n) -Bz sparse(n,n); Bz  sparse(n,n) sparse(n,n); -By Bx sparse(n,n)];    
elseif strcmp(modetype, 'TE')
    F = [sparse(n,n) sparse(n,n) sparse(n,n); Fz sparse(n,n) sparse(n,n); -Fy Fx sparse(n,n)];
    B = [sparse(n,n) -Bz By; Bz  sparse(n,n) -Bx; sparse(n,n) sparse(n,n) sparse(n,n)];
end

double_rotor = [sparse(3*n, 3*n) F/1i; B/1i sparse(3*n, 3*n)];

%FORWARD AND DIAG MATRIX
function Dout = D(n)
    Dout = speye(n);
end
function Dpout = Dp(n)
    Dpout = sparse(1:n-1,2:n,ones(1,n-1),n,n);
end
function fiout = fiX(n,k)
    fiout = Dp(n) - D(n);
    fiout(n, 1) = exp(1i*k(1));
end
function fiout = fiY(n,k)
    fiout = Dp(n) - D(n);
    fiout(n, 1) = exp(1i*k(2));
end
function fiout = fiZ(n,k)
    fiout = Dp(n) - D(n);
    fiout(n, 1) = exp(1i*k(3));
end
end

