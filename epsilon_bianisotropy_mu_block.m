function eps_chi_chi_mu = epsilon_bianisotropy_mu_block(epstype, biantype, mutype, lambda2)
%THIS FUNCTION CALCULATES THE EPSILON MATRIX ELEMENT THROUGHT THE ENTIRE
%DOMAIN - 

global nx ny nz u1 u2 u3 NO DN NOUT file R c gain Delta modetype;

%calculates transition matrix P: standard basis -> basis of new basis
%vectors u1, u2, u3 and it`s inverse
Pmatrix = [u1', u2', u3'];
PmatrixINV = inv(Pmatrix);

n = nx*ny*nz;

plotdir(epstype);

Di = zeros(n,1);
Dj = zeros(n,1);

Exxv = zeros(n,1); Exyv = zeros(n,1); Exzv = zeros(n,1); 
Eyxv = zeros(n,1); Eyyv = zeros(n,1); Eyzv = zeros(n,1);
Ezxv = zeros(n,1); Ezyv = zeros(n,1);  Ezzv = zeros(n,1);

Bxxv = zeros(n,1); Bxyv = zeros(n,1); Bxzv = zeros(n,1); 
Byxv = zeros(n,1); Byyv = zeros(n,1); Byzv = zeros(n,1);
Bzxv = zeros(n,1); Bzyv = zeros(n,1); Bzzv = zeros(n,1);

Mxxv = zeros(n,1); Mxyv = zeros(n,1); Mxzv = zeros(n,1); 
Myxv = zeros(n,1); Myyv = zeros(n,1); Myzv = zeros(n,1);
Mzxv = zeros(n,1); Mzyv = zeros(n,1); Mzzv = zeros(n,1);

Rxpi = zeros(2*n,1); Rypi = zeros(2*n,1); Rzpi = zeros(2*n,1);
Rxpj = zeros(2*n,1); Rypj = zeros(2*n,1); Rzpj = zeros(2*n,1);
Rxpv = zeros(2*n,1); Rypv = zeros(2*n,1); Rzpv = zeros(2*n,1);


%FILL THE ARRAYS OF EPSILON COMPONENTS

for i = 1:nx
    for j = 1:ny
        for k = 1:nz      
            D1 = eblock(epstype,i,j,k); %epsilon
            if strcmp(mutype, 'IDENTITY')
                D2 = speye(3);
            else
                D2 = eblock(mutype,i,j,k); %mu
            end
            D3 = bianblock(biantype,i,j,k); %bianisotropy
            
            % epsilon and mu blocks are all diagonal matrices
            Di(ind(i,j,k)) = ind(i,j,k);
            Dj(ind(i,j,k)) = ind(i,j,k);
            
            %dielectric tensor
            Exxv(ind(i,j,k)) = D1(1,1); Mxxv(ind(i,j,k)) = D2(1,1); Bxxv(ind(i,j,k)) = D3(1,1); 
            Exyv(ind(i,j,k)) = D1(1,2); Mxyv(ind(i,j,k)) = D2(1,2); Bxyv(ind(i,j,k)) = D3(1,2);
            Exzv(ind(i,j,k)) = D1(1,3); Mxzv(ind(i,j,k)) = D2(1,3); Bxzv(ind(i,j,k)) = D3(1,3);
                
            Eyxv(ind(i,j,k)) = D1(2,1); Myxv(ind(i,j,k)) = D2(2,1); Byxv(ind(i,j,k)) = D3(2,1);
            Eyyv(ind(i,j,k)) = D1(2,2); Myyv(ind(i,j,k)) = D2(2,2); Byyv(ind(i,j,k)) = D3(2,2); 
            Eyzv(ind(i,j,k)) = D1(2,3); Myzv(ind(i,j,k)) = D2(2,3); Byzv(ind(i,j,k)) = D3(2,3);
                
            Ezxv(ind(i,j,k)) = D1(3,1); Mzxv(ind(i,j,k)) = D2(3,1); Bzxv(ind(i,j,k)) = D3(3,1);
            Ezyv(ind(i,j,k)) = D1(3,2); Mzyv(ind(i,j,k)) = D2(3,2); Bzyv(ind(i,j,k)) = D3(3,2);
            Ezzv(ind(i,j,k)) = D1(3,3); Mzzv(ind(i,j,k)) = D2(3,3); Bzzv(ind(i,j,k)) = D3(3,3); 
        end
    end
end

%interpolation matrices: shift a field to the nearby points, see Raymond, Finite-Difference Frequency-Domain Algorithm for Modeling Electromagnetic Scattering from General Anisotropic Objects
for i = 1:nx
    for j = 1:ny
        for k = 1:nz  
            Rxpi(2*ind(i,j,k)-1) = ind(i,j,k); %diagonalni
            Rxpj(2*ind(i,j,k)-1) = ind(i,j,k); %diagonalni
            Rypi(2*ind(i,j,k)-1) = ind(i,j,k); %diagonalni
            Rypj(2*ind(i,j,k)-1) = ind(i,j,k); %diagonalni
            Rzpi(2*ind(i,j,k)-1) = ind(i,j,k); %diagonalni
            Rzpj(2*ind(i,j,k)-1) = ind(i,j,k); %diagonalni
            
            Rxpi(2*ind(i,j,k)-0) = ind(i,j,k); %izvendiagonalni
            Rxpj(2*ind(i,j,k)-0) = ind(i+1,j,k); %izvendiagonalni
            Rypi(2*ind(i,j,k)-0) = ind(i,j,k); %izvendiagonalni
            Rypj(2*ind(i,j,k)-0) = ind(i,j+1,k); %izvendiagonalni
            Rzpi(2*ind(i,j,k)-0) = ind(i,j,k); %izvendiagonalni
            Rzpj(2*ind(i,j,k)-0) = ind(i,j,k+1); %izvendiagonalni
                       
            if i == nx
                Rxpv(2*ind(i,j,k)-1) = 1; %maybe2??;       %diagonalni
                Rxpv(2*ind(i,j,k)-0) = 0;       %diagonalni
            else
                Rxpv(2*ind(i,j,k)-1) = 1;       %diagonalni
                Rxpv(2*ind(i,j,k)-0) = 1;       %izvendiagonalni
            end
            
            if j == ny
                Rypv(2*ind(i,j,k)-1) = 1; %maybe2??;       %diagonalni
                Rypv(2*ind(i,j,k)-0) = 0;       %diagonalni
            else
                Rypv(2*ind(i,j,k)-1) = 1;       %diagonalni
                Rypv(2*ind(i,j,k)-0) = 1;       %izvendiagonalni
            end
            
            if k == nz
                Rzpv(2*ind(i,j,k)-1) = 1; %maybe2??;       %diagonalni
                Rzpv(2*ind(i,j,k)-0) = 0;       %diagonalni
            else
                Rzpv(2*ind(i,j,k)-1) = 1;       %diagonalni
                Rzpv(2*ind(i,j,k)-0) = 1;       %izvendiagonalni
            end
            
        end
    end
end

Exx = sparse(Di,Dj,Exxv);
Exy = sparse(Di,Dj,Exyv);
Exz = sparse(Di,Dj,Exzv);

Eyx = sparse(Di,Dj,Eyxv);
Eyy = sparse(Di,Dj,Eyyv);
Eyz = sparse(Di,Dj,Eyzv);

Ezx = sparse(Di,Dj,Ezxv);
Ezy = sparse(Di,Dj,Ezyv);
Ezz = sparse(Di,Dj,Ezzv);

Mxx = sparse(Di,Dj,Mxxv);
Mxy = sparse(Di,Dj,Mxyv);
Mxz = sparse(Di,Dj,Mxzv);

Myx = sparse(Di,Dj,Myxv);
Myy = sparse(Di,Dj,Myyv);
Myz = sparse(Di,Dj,Myzv);

Mzx = sparse(Di,Dj,Mzxv);
Mzy = sparse(Di,Dj,Mzyv);
Mzz = sparse(Di,Dj,Mzzv);

Bxx = sparse(Di,Dj,Bxxv);
Bxy = sparse(Di,Dj,Bxyv);
Bxz = sparse(Di,Dj,Bxzv);

Byx = sparse(Di,Dj,Byxv);
Byy = sparse(Di,Dj,Byyv);
Byz = sparse(Di,Dj,Byzv);

Bzx = sparse(Di,Dj,Bzxv);
Bzy = sparse(Di,Dj,Bzyv);
Bzz = sparse(Di,Dj,Bzzv);


Rxp = 1/2*sparse(Rxpi,Rxpj,Rxpv);
Ryp = 1/2*sparse(Rypi,Rypj,Rypv);
Rzp = 1/2*sparse(Rzpi,Rzpj,Rzpv);

Rxp = Rxp(1:nx*ny*nz,1:nx*ny*nz); 
Ryp = Ryp(1:nx*ny*nz,1:nx*ny*nz); 
Rzp = Rzp(1:nx*ny*nz,1:nx*ny*nz); 

Rxm = Rxp';
Rym = Ryp';
Rzm = Rzp';

mat0 = sparse(n,n);
if strcmp(modetype, 'ALL') 
    Eout = [        Exx, Rxp*Rym*Exy,     Rxp*Rzm*Exz; 
            Ryp*Rxm*Eyx,         Eyy,     Ryp*Rzm*Eyz; 
            Rzp*Rxm*Ezx, Rzp*Rym*Ezy,             Ezz];


    Mout = [        Mxx, Rxm*Ryp*Mxy,     Rxm*Rzp*Mxz; 
            Rym*Rxp*Myx,         Myy,     Rym*Rzp*Myz;
            Rzm*Rxp*Mzx, Rzm*Ryp*Mzy,             Mzz];

    BoutE = conj([Rxm*Ryp*Rzp*Bxx, Rxp*Rym*Rzp*Byx, Rxp*Ryp*Rzm*Bzx;
                  Rxm*Ryp*Rzp*Bxy, Rxp*Rym*Rzp*Byy, Rxp*Ryp*Rzm*Bzy;
                  Rxm*Ryp*Rzp*Bxz, Rxp*Rym*Rzp*Byz, Rxp*Ryp*Rzm*Bzz]);

    BoutH = [Rxp*Rym*Rzm*Bxx, Rxm*Ryp*Rzm*Bxy, Rxm*Rym*Rzp*Bxz;
             Rxp*Rym*Rzm*Byx, Rxm*Ryp*Rzm*Byy, Rxm*Rym*Rzp*Byz;
             Rxp*Rym*Rzm*Bzx, Rxm*Ryp*Rzm*Bzy, Rxm*Rym*Rzp*Bzz];
elseif strcmp(modetype, 'TM') 
    Eout = [mat0, mat0,     Rxp*Rzm*Exz; 
            mat0, mat0,     Ryp*Rzm*Eyz; 
            mat0, mat0,             Ezz];


    Mout = [mat0,        mat0,         mat0; 
            mat0,        mat0,         mat0;
            Rzm*Rxp*Mzx, Rzm*Ryp*Mzy,  Mzz];

    BoutE = conj([mat0, mat0, Rxp*Ryp*Rzm*Bzx;
                  mat0, mat0, Rxp*Ryp*Rzm*Bzy;
                  mat0, mat0, Rxp*Ryp*Rzm*Bzz]);

    BoutH = [mat0,            mat0,            mat0;
             mat0,            mat0,            mat0;
             Rxp*Rym*Rzm*Bzx, Rxm*Ryp*Rzm*Bzy, Rxm*Rym*Rzp*Bzz];
elseif strcmp(modetype, 'TE') 
    Eout = [mat0,        mat0,         mat0; 
            mat0,        mat0,         mat0;
            Rzp*Rxm*Ezx, Rzp*Rym*Ezy,  Ezz];


    Mout = [mat0, mat0,     Rxm*Rzp*Mxz; 
            mat0, mat0,     Rym*Rzp*Myz;
            mat0, mat0,             Mzz];

    BoutE = conj([mat0,            mat0,            mat0;
                  mat0,            mat0,            mat0;
                  Rxm*Ryp*Rzp*Bxz, Rxp*Rym*Rzp*Byz, Rxp*Ryp*Rzm*Bzz]);

    BoutH = [mat0, mat0, Rxm*Rym*Rzp*Bxz;
             mat0, mat0, Rxm*Rym*Rzp*Byz;
             mat0, mat0, Rxm*Rym*Rzp*Bzz];
end
     
eps_chi_chi_mu = [-Eout BoutH; BoutE Mout];
     
%define space dependance of refractive index
function f = nord(i,j,k)
    TW = tukeywin(200, 0.4);
    if r(i,j,nz/2,c) < R
        f = NO; %NOUT + (NO - NOUT)*TW(round(100*(1-r(i,j,k,c)/R)+1));
    else
        f = NOUT;
    end
end
function f = dn(i,j,k)
  f = DN;
end

function f = ind(i,j,k)
    f = (k-1)*nx*ny + (j-1)*nx + (i-1) + 1;
end

%SELECT DIRECTOR FIELD

function [f1,f2,f3] = dir(dir_type, ri0,rj0,rk0)
    %transform to standard coordinates
    standard_index_vector = PmatrixINV * [ri0;rj0;rk0];
    ri = standard_index_vector(1); rj = standard_index_vector(2); rk = standard_index_vector(3);
    
%RANDOM
    if strcmp(dir_type, 'RANDOM')
        f11 = rand;
        f21 = rand;
        f31 = rand;
        f1 = f11/sqrt(f11^2+f21^2+f31^2);
        f2 = f21/sqrt(f11^2+f21^2+f31^2);
        f3 = f31/sqrt(f11^2+f21^2+f31^2);
%ISOTROPIC
    elseif strcmp(dir_type, 'ISOTROPIC')
           f1 = 1;
           f2 = 0;
           f3 = 0; 
%RADIAL DROPLET
    elseif strcmp(dir_type, 'ISOSPHERE')
       if r(ri,rj,rk,c) < R 
           f1 = 1;
           f2 = 0;
           f3 = 0;           
       else 
           f1 = 0;
           f2 = 0;
           f3 = 0;
       end        
%ZERO 
    elseif strcmp(dir_type, 'ZERO')
        f1 = 0;
        f2 = 0;
        f3 = 0;
% HELICONICAL
    elseif strcmp(dir_type, 'HELICONIC')
        f1 = cos(pi/4);
        f2 = sin(2*pi*ri/(nx))*sin(pi/4);
        f3 = cos(2*pi*ri/(nx))*sin(pi/4);
%HELICONICAL XY DIAGONAL
    elseif strcmp(dir_type,'HELICONICXY')
        pitch = nx*sqrt(2)/2;
        f1 = cos(pi/4)*cos(pi/4) - sin(phi(ri*cos(pi/4)+rj*(sin(pi/4)),pitch)) * sin(pi/4) * sin(pi/4);
        f2 = cos(pi/4)*sin(pi/4) + sin(phi(ri*cos(pi/4)+rj*(sin(pi/4)),pitch)) * sin(pi/4) * cos(pi/4);
        f3 = cos(phi(ri*cos(pi/4)+rj*(sin(pi/4)),pitch))*sin(pi/4);
%HELICONICAL XZ DIAGONAL
    elseif strcmp(dir_type, 'HELICONICXZ')
        theta = pi/4;
        pitch = nx*sqrt(2)/2;
        f1 = cos(theta)*cos(pi/4) - cos(phi(ri*cos(pi/4)+rk*(sin(pi/4)),pitch))*sin(theta) * sin(pi/4);
        f2 = sin(phi(ri*cos(pi/4)+rk*(sin(pi/4)),pitch)) * sin(theta);
        f3 = cos(theta)*sin(pi/4) + cos(phi(ri*cos(pi/4)+rk*(sin(pi/4)),pitch))*sin(theta) * cos(pi/4);
%HELICONICAL YZ DIAGONAL
    elseif strcmp(dir_type, 'HELICONICYZ')
        theta = pi/4;
        pitch = ny*sqrt(2)/2;
        f1 = cos(theta)*cos(pi/4) - cos(phi(rj*cos(pi/4)+rk*(sin(pi/4)),pitch))*sin(theta) * sin(pi/4);
        f2 = sin(phi(rj*cos(pi/4)+rk*(sin(pi/4)),pitch)) * sin(theta);
        f3 = cos(theta)*sin(pi/4) + cos(phi(rj*cos(pi/4)+rk*(sin(pi/4)),pitch))*sin(theta) * cos(pi/4);        
%RADIAL DROPLET
    elseif strcmp(dir_type, 'RADIALD')
       if r(ri,rj,rk,c) < R 
           f1 = -(ri-c(1))/(r(ri,rj,rk,c)+1E-12);
           f2 = -(rj-c(2))/(r(ri,rj,rk,c)+1E-12);
           f3 = -(rk-c(3))/(r(ri,rj,rk,c)+1E-12);           
       else 
           f1 = 0;
           f2 = 0;
           f3 = 0;
       end
%QUARTER OF A DROPLET
    elseif strcmp(dir_type, 'RADIALD_QTR')
       if r_qrt(ri,rj,rk) < R 
           f1 = -(ri)/(r_qrt(ri,rj,rk)+1E-12);
           f2 = -(rj)/(r_qrt(ri,rj,rk)+1E-12);
           f3 = -(rk)/(r_qrt(ri,rj,rk)+1E-12);           
       else 
           f1 = 0;
           f2 = 0;
           f3 = 0;
       end
%ESCAPED CYLINDER GUIDE
    elseif strcmp(dir_type, 'ESCAPEDC')
       if sqrt((ri-nx/2)^2 + (rj-ny/2)^2) < nx/4 
           f1 = -sin(2*pi*(ri-nx/2)/nx);
           f2 = -sin(2*pi*(rj-ny/2)/nx);
           f3 = sqrt(1-f1^2-f2^2);           
       else 
           f1 = -(ri-nx/2)/sqrt((ri-nx/2)^2 + (rj-ny/2)^2 +1E-12);
           f2 = -(rj-ny/2)/sqrt((ri-nx/2)^2 + (rj-ny/2)^2 +1E-12);
           f3 = 0;
       end
%BIPOLAR DROPLET
    elseif strcmp(dir_type, 'BIPOLAR')
       if r(ri,rj,rk,c) < R 
           th = (rk-c(3))/(R*r(ri/R,rj/R,rk/R,c/R));
           ph = atan2((rj-c(2))/R,(ri-c(1))/R);
           f1 = cos(th)*cos(ph) / sqrt(cos(th)^2 + sin(th)^2);
           f2 = cos(th)*sin(ph) / (cos(th)^2 + sin(th)^2);
           f3 = -sin(th) / sqrt(cos(th)^2 + sin(th)^2);       
       else 
           f1 = 0;
           f2 = 0;
           f3 = 0;
       end
%BIPOLAR DROPLET in x direction
     elseif strcmp(dir_type, 'BIPOLARX')
       if r(ri,rj,rk,c) < R 
           th = (ri-c(3))/(R*r(ri/R,rj/R,rk/R,c/R));
           ph = atan2((rj-c(2))/R,(rk-c(3))/R);
           f3 = -cos(th)*cos(ph) / sqrt(cos(th)^2 + sin(th)^2);
           f2 = cos(th)*sin(ph) / (cos(th)^2 + sin(th)^2);
           f1 = -sin(th) / sqrt(cos(th)^2 + sin(th)^2);       
       else 
           f1 = 0;
           f2 = 0;
           f3 = 0;
       end
%CYLINDER
       elseif strcmp(dir_type, 'CYLINDER')
       if sqrt((ri-c(1))^2 + (rj-c(2))^2) < R 
           f1 = 0;
           f2 = 1;
           f3 = 0;           
       else 
           f1 = 0;
           f2 = 0;
           f3 = 0;
       end
%read DIR from file   
     elseif strcmp(dir_type, 'FILE')  
         f1 = dirarray(3*ind(ri,rj,rk)-2);
         f2 = dirarray(3*ind(ri,rj,rk)-1);
         f3 = dirarray(3*ind(ri,rj,rk)-0);
    else
        disp('ERROR: Unknown analytical director field');
        return;
    end    
end  


%phi angle for heliconics and cholesterics
function f = phi(l,pitch)
    f = 2*pi*l/(pitch);
end

%radius for droplets
function f = r(i,j,k,c)
       f = sqrt((i-c(1))^2 + (j-c(2))^2 + (k-c(3))^2);
end

%radius for droplets
function f = r_qrt(i,j,k)
       f = sqrt((i)^2 + (j)^2 + (k)^2);
end

%calculates eps_o and deps
function f = epsO(i,j,k)
    if abs(gain) > 0 && r(i,j,k,c) < R
        f = ((nord(i,j,k)+dn(i,j,k))^2 + 2*nord(i,j,k)^2)/3 - gain^2/4;
    else
        f = ((nord(i,j,k)+dn(i,j,k))^2 + 2*nord(i,j,k)^2)/3;
    end
end
function f = deps(i,j,k)
    f = (nord(i,j,k)+dn(i,j,k))^2 - nord(i,j,k)^2;
%     TW = tukeywin(2*nx/4, 0.4);
%     if r(i,j,k)<nx/4
%         f = (nord(i,j,k)+TW(round(r(i,j,k))+nx/4)*dn(i,j,k))^2 - nord(i,j,k)^2;
%     else
%         f = (nord(i,j,k)+dn(i,j,k))^2 - nord(i,j,k)^2;
%     end
end


%calculates components of epsilon from director field, accounts for Yee
%lattice shifts
function eps = eblock(dir_type, i,j,k)    
    %i+1/2 for exx eyx ezx   
    [dir1,dir2,dir3]=dir(dir_type, i +1/2,j,k);

    Qxx = (dir1*dir1 * 3 -1)/2;
    Qxy = (dir1*dir2 * 3)/2;
    Qxz = (dir1*dir3 * 3)/2;
    Qyy = (dir2*dir2 * 3 -1)/2;
    Qyz = (dir2*dir3 * 3)/2;
    Qzz = (dir3*dir3 * 3 -1)/2;
    
    Qten = [Qxx Qxy Qxz; Qxy Qyy Qyz; Qxz Qyz Qzz];
    
    epsX = (deps(i +1/2,j,k) * Qten + (epsO(i +1/2,j,k)+deps(i +1/2,j,k)/3) * eye (3)); 

    %j+1/2 for exy eyy ezy
    [dir1,dir2,dir3]=dir(dir_type, i,j +1/2,k);
    Qxx = (dir1*dir1 * 3 -1)/2;
    Qxy = (dir1*dir2 * 3)/2;
    Qxz = (dir1*dir3 * 3)/2;
    Qyy = (dir2*dir2 * 3 -1)/2;
    Qyz = (dir2*dir3 * 3)/2;
    Qzz = (dir3*dir3 * 3 -1)/2;
    
    Qten = [Qxx Qxy Qxz; Qxy Qyy Qyz; Qxz Qyz Qzz];
    
    epsY = (deps(i,j +1/2,k) * Qten + (epsO(i,j +1/2,k)+deps(i,j +1/2,k)/3) * eye (3)); 
   
    
    %k+1/2 for exz eyz ezz
    [dir1,dir2,dir3]=dir(dir_type, i,j,k+1/2);
    Qxx = (dir1*dir1 * 3 -1)/2;
    Qxy = (dir1*dir2 * 3)/2;
    Qxz = (dir1*dir3 * 3)/2;
    Qyy = (dir2*dir2 * 3 -1)/2;
    Qyz = (dir2*dir3 * 3)/2;
    Qzz = (dir3*dir3 * 3 -1)/2;
    
    Qten = [Qxx Qxy Qxz; Qxy Qyy Qyz; Qxz Qyz Qzz];
    
    epsZ = (deps(i,j,k +1/2) * Qten + (epsO(i,j,k +1/2)+deps(i,j,k +1/2)/3) * eye (3));
    
    eps = [epsX(1,1) epsY(1,2) epsZ(1,3); epsX(2,1) epsY(2,2) epsZ(2,3); epsX(3,1) epsY(3,2) epsZ(3,3)];
    eps = Pmatrix*eps*PmatrixINV;
end

function bian = bianblock(bian_type, ri0,rj0,rk0)
    standard_index_vector = PmatrixINV * [ri0;rj0;rk0];
    ri = standard_index_vector(1); rj = standard_index_vector(2); rk = standard_index_vector(3);
    
%ZERO
    if strcmp(bian_type, 'ZERO')
        bian = sparse(3, 3);
%ISOTROPIC
    elseif strcmp(bian_type, 'KHANIKAEV')
        if r(ri,rj,rk,c) < R 
            bian = [0 1i*Delta 0; -1i*Delta 0 0; 0 0 0];
        else
            bian = sparse(3, 3);
        end
    end
    bian = Pmatrix * bian * PmatrixINV;
end

function [] = plotdir(dir_type)
    u1=zeros(nx,ny);
    v1=zeros(nx,ny);
    u2=zeros(nx,nz);
    v2=zeros(nx,nz);
    u3=zeros(ny,nz);
    v3=zeros(ny,nz);
    nordplot = zeros(nx);
    dnplot = zeros(nx);
    nordplot2 = zeros(nx);
    
     for ii = 1:nx 
        for jj = 1:ny
            for kk = 1:nz
%                 [u1(i,j),v1(i,j),~] = dir(i,j,nz-1);
%                 [u2(i,k),~,v2(i,k)] = dir(i,1,k);
%                 [~,u3(j,k),v3(j,k)] = dir(nx - 1,j,k); 
                [u1(ii,jj),v1(ii,jj),~] = dir(dir_type,ii,jj,(nz+1)/2);
                [u2(ii,kk),~,v2(ii,kk)] = dir(dir_type,ii,(nz+1)/2+1,kk);
                [~,u3(jj,kk),v3(jj,kk)] = dir(dir_type,(nz+1)/2+1,jj,kk); 

            end
        end
     end

    [x1,y1]=meshgrid(1:nx,1:ny);
    [x2,z2]=meshgrid(1:nx,1:nz);
    [y3,z3]=meshgrid(1:ny,1:nz);
    
    for ii=1:nx
        nordplot(ii) = nord(ii,3*ny/4,3*nz/4);
        dnplot(ii) = dn(ii,3*ny/4,3*nz/4);
        nordplot2(ii) =nord(ii,1,1);
    end

    f1 = figure;
        set(f1,'visible','off');
        quiver(x1,y1,u1',v1');
        pbaspect([1 1 1]);
    f2 = figure;
        set(f2,'visible','off');
        quiver(x2,z2,u2',v2');
        pbaspect([1 1 1]);
    f3 = figure;
        set(f3,'visible','off');
        quiver(y3,z3,u3',v3');
        pbaspect([1 1 1]);
    f4 = figure;
        set(f4,'visible','off');    
        plot(nordplot);
    f5 = figure;
        set(f5,'visible','off');   
        plot(dnplot);   
    f6 = figure;
        set(f6,'visible','off');   
        plot(nordplot2);
        
     out = strcat(file,'/dir');
     print(f1,strcat(out,'XY'),'-dpng');
     print(f2,strcat(out,'XZ'),'-dpng');
     print(f3,strcat(out,'YZ'),'-dpng');
     print(f4,strcat(out,'nord'),'-dpng');
     print(f5,strcat(out,'dn'),'-dpng');
     print(f6,strcat(out,'nord2'),'-dpng');
end

end