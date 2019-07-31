function [Eout,Mout] = epsilon_S_block(pdir,lambda2)
%THIS FUNCTION CALCULATES THE EPSILON MATRIX ELEMENT THROUGHT THE ENTIRE
%DOMAIN - 

global nx ny nz NO DN NOUT dirfield file R c DPMLs DPMLe gain;

%PLOT DIRECTOR FIELD
if strcmp(pdir, 'YES')
    
    if strcmp(dirfield, 'FILE')  
        dirarray = read_dir('/data/Matlab/data/141128_22_input_M3_res_dir.raw');
    end
    
    u1=zeros(nx,ny);
    v1=zeros(nx,ny);
    u2=zeros(nx,nz);
    v2=zeros(nx,nz);
    u3=zeros(ny,nz);
    v3=zeros(ny,nz);
    nordplot = zeros(nx);
    dnplot = zeros(nx);
    nordplot2 = zeros(nx);
    
    exportdir(dirfield);
    
     for i = 1:nx 
        for j = 1:ny
            for k = 1:nz
%                 [u1(i,j),v1(i,j),~] = dir(i,j,nz-1);
%                 [u2(i,k),~,v2(i,k)] = dir(i,1,k);
%                 [~,u3(j,k),v3(j,k)] = dir(nx - 1,j,k); 
                [u1(i,j),v1(i,j),~] = dir(i,j,(nz+1)/2);
                [u2(i,k),~,v2(i,k)] = dir(i,(nz+1)/2+1,k);
                [~,u3(j,k),v3(j,k)] = dir((nz+1)/2+1,j,k); 

            end
        end
     end

    [x1,y1]=meshgrid(1:nx,1:ny);
    [x2,z2]=meshgrid(1:nx,1:nz);
    [y3,z3]=meshgrid(1:ny,1:nz);
    
    for i=1:nx
        nordplot(i) = nord(i,3*ny/4,3*nz/4);
        dnplot(i) = dn(i,3*ny/4,3*nz/4);
        nordplot2(i) =nord(i,1,1);
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
        
    if isempty(file) == 0
        out = strcat(file,'/dir');
        print(f1,strcat(out,'XY'),'-dpng');
        print(f2,strcat(out,'XZ'),'-dpng');
        print(f3,strcat(out,'YZ'),'-dpng');
        print(f4,strcat(out,'nord'),'-dpng');
        print(f5,strcat(out,'dn'),'-dpng');
        print(f6,strcat(out,'nord2'),'-dpng');
   end
end

n = nx*ny*nz;

Di = zeros(n,1);
Dj = zeros(n,1);

Exxv = zeros(n,1);
Exyv = zeros(n,1); 
Exzv = zeros(n,1); 

Eyxv = zeros(n,1);
Eyyv = zeros(n,1); 
Eyzv = zeros(n,1);

Ezxv = zeros(n,1);
Ezyv = zeros(n,1); 
Ezzv = zeros(n,1);

Mxxv = zeros(n,1);
Myyv = zeros(n,1); 
Mzzv = zeros(n,1);

Rxpi = zeros(2*n,1);
Rypi = zeros(2*n,1);
Rzpi = zeros(2*n,1);
Rxpj = zeros(2*n,1);
Rypj = zeros(2*n,1);
Rzpj = zeros(2*n,1);
Rxpv = zeros(2*n,1);
Rypv = zeros(2*n,1);
Rzpv = zeros(2*n,1);


%FILL THE ARRAYS OF EPSILON COMPONENTS
for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            
            [D1,DM1] = eblock(i,j,k);
            
            % epsilon and mu blocks are all diagonal matrices
            Di(ind(i,j,k)) = ind(i,j,k);
            Dj(ind(i,j,k)) = ind(i,j,k);
            
            %dielectric tensor
            Exxv(ind(i,j,k)) = D1(1,1); 
            Exyv(ind(i,j,k)) = D1(1,2);
            Exzv(ind(i,j,k)) = D1(1,3);
                
            Eyxv(ind(i,j,k)) = D1(2,1);
            Eyyv(ind(i,j,k)) = D1(2,2); 
            Eyzv(ind(i,j,k)) = D1(2,3);
                
            Ezxv(ind(i,j,k)) = D1(3,1);
            Ezyv(ind(i,j,k)) = D1(3,2);
            Ezzv(ind(i,j,k)) = D1(3,3); 
              
            %permeability matrix with PML
            Mxxv(ind(i,j,k)) = DM1(1,1); 
            Myyv(ind(i,j,k)) = DM1(2,2); 
            Mzzv(ind(i,j,k)) = DM1(3,3); 
            
            %espilon interpolation matrices
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
Myy = sparse(Di,Dj,Myyv);
Mzz = sparse(Di,Dj,Mzzv);

Rxp = 1/2*sparse(Rxpi,Rxpj,Rxpv);
Ryp = 1/2*sparse(Rypi,Rypj,Rypv);
Rzp = 1/2*sparse(Rzpi,Rzpj,Rzpv);

Rxp = Rxp(1:nx*ny*nz,1:nx*ny*nz); 
Ryp = Ryp(1:nx*ny*nz,1:nx*ny*nz); 
Rzp = Rzp(1:nx*ny*nz,1:nx*ny*nz); 

Rxm = Rxp';
Rym = Ryp';
Rzm = Rzp';

Eout = [Exx Rxp*Rym*Exy Rxp*Rzm*Exz; Ryp*Rxm*Eyx Eyy Ryp*Rzm*Eyz; Rzp*Rxm*Ezx Rzp*Rym*Ezy Ezz];
Mout = [Mxx sparse(n,n) sparse(n,n); sparse(n,n) Myy sparse(n,n); sparse(n,n) sparse(n,n) Mzz];

%define space dependance of refractive index
function f = nord(i,j,k)
    TW = tukeywin(200, 0.4);
    if r(i,j,k,c) < R
        f = NO; %NOUT + (NO - NOUT)*TW(round(100*(1-r(i,j,k,c)/R)+1));
    else
        f = NOUT;
    end
end
function f = dn(i,j,k)
  f = DN;
%     TW = tukeywin(200, 0.4);
%     if r(i,j,k,c) < R
%         f = DN; %*TW(round(100*(1-r(i,j,k,c)/R)+1));
%     else
%         f = 0;
%     end
end

function f = ind(i,j,k)
    f = (k-1)*nx*ny + (j-1)*nx + (i-1) + 1;
end

%SELECT DIRECTOR FIELD

function [f1,f2,f3] = dir(i,j,k)
%RANDOM
    if strcmp(dirfield, 'RANDOM')
        f11 = rand;
        f21 = rand;
        f31 = rand;
        f1 = f11/sqrt(f11^2+f21^2+f31^2);
        f2 = f21/sqrt(f11^2+f21^2+f31^2);
        f3 = f31/sqrt(f11^2+f21^2+f31^2);
%ISOTROPIC
    elseif strcmp(dirfield, 'ISOTROPIC')
           f1 = 1;
           f2 = 1;
           f3 = 1; 
%            f1 = -(i-c(1))/(sqrt((i-c(1))^2 + (j-c(2))^2 )+1E-12);
%            f2 = -(j-c(2))/(sqrt((i-c(1))^2 + (j-c(2))^2 )+1E-12);
%            f3 = 0; 
%RADIAL DROPLET
    elseif strcmp(dirfield, 'ISOSPHERE')
       if r(i,j,k,c) < R 
           f1 = 1;
           f2 = 0;
           f3 = 0;           
       else 
           f1 = 0;
           f2 = 0;
           f3 = 0;
       end        
%ZERO 
    elseif strcmp(dirfield, 'ZERO')
        f1 = 0;
        f2 = 0;
        f3 = 0;
% HELICONICAL
    elseif strcmp(dirfield, 'HELICONIC')
        f1 = cos(pi/4);
        f2 = sin(2*pi*i/(nx))*sin(pi/4);
        f3 = cos(2*pi*i/(nx))*sin(pi/4);
%HELICONICAL XY DIAGONAL
    elseif strcmp(dirfield,'HELICONICXY')
        pitch = nx*sqrt(2)/2;
        f1 = cos(pi/4)*cos(pi/4) - sin(phi(i*cos(pi/4)+j*(sin(pi/4)),pitch)) * sin(pi/4) * sin(pi/4);
        f2 = cos(pi/4)*sin(pi/4) + sin(phi(i*cos(pi/4)+j*(sin(pi/4)),pitch)) * sin(pi/4) * cos(pi/4);
        f3 = cos(phi(i*cos(pi/4)+j*(sin(pi/4)),pitch))*sin(pi/4);
%HELICONICAL XZ DIAGONAL
    elseif strcmp(dirfield, 'HELICONICXZ')
        theta = pi/4;
        pitch = nx*sqrt(2)/2;
        f1 = cos(theta)*cos(pi/4) - cos(phi(i*cos(pi/4)+k*(sin(pi/4)),pitch))*sin(theta) * sin(pi/4);
        f2 = sin(phi(i*cos(pi/4)+k*(sin(pi/4)),pitch)) * sin(theta);
        f3 = cos(theta)*sin(pi/4) + cos(phi(i*cos(pi/4)+k*(sin(pi/4)),pitch))*sin(theta) * cos(pi/4);
%HELICONICAL YZ DIAGONAL
    elseif strcmp(dirfield, 'HELICONICYZ')
        theta = pi/4;
        pitch = ny*sqrt(2)/2;
        f1 = cos(theta)*cos(pi/4) - cos(phi(j*cos(pi/4)+k*(sin(pi/4)),pitch))*sin(theta) * sin(pi/4);
        f2 = sin(phi(j*cos(pi/4)+k*(sin(pi/4)),pitch)) * sin(theta);
        f3 = cos(theta)*sin(pi/4) + cos(phi(j*cos(pi/4)+k*(sin(pi/4)),pitch))*sin(theta) * cos(pi/4);        
%RADIAL DROPLET
    elseif strcmp(dirfield, 'RADIALD')
       if r(i,j,k,c) < R 
           f1 = -(i-c(1))/(r(i,j,k,c)+1E-12);
           f2 = -(j-c(2))/(r(i,j,k,c)+1E-12);
           f3 = -(k-c(3))/(r(i,j,k,c)+1E-12);           
       else 
           f1 = 0;
           f2 = 0;
           f3 = 0;
       end
%QUARTER OF A DROPLET
    elseif strcmp(dirfield, 'RADIALD_QTR')
       if r_qrt(i,j,k) < R 
           f1 = -(i)/(r_qrt(i,j,k)+1E-12);
           f2 = -(j)/(r_qrt(i,j,k)+1E-12);
           f3 = -(k)/(r_qrt(i,j,k)+1E-12);           
       else 
           f1 = 0;
           f2 = 0;
           f3 = 0;
       end
%ESCAPED CYLINDER GUIDE
    elseif strcmp(dirfield, 'ESCAPEDC')
       if sqrt((i-nx/2)^2 + (j-ny/2)^2) < nx/4 
           f1 = -sin(2*pi*(i-nx/2)/nx);
           f2 = -sin(2*pi*(j-ny/2)/nx);
           f3 = sqrt(1-f1^2-f2^2);           
       else 
           f1 = -(i-nx/2)/sqrt((i-nx/2)^2 + (j-ny/2)^2 +1E-12);
           f2 = -(j-ny/2)/sqrt((i-nx/2)^2 + (j-ny/2)^2 +1E-12);
           f3 = 0;
       end
%BIPOLAR DROPLET
    elseif strcmp(dirfield, 'BIPOLAR')
       if r(i,j,k,c) < R 
           th = (k-c(3))/(R*r(i/R,j/R,k/R,c/R));
           ph = atan2((j-c(2))/R,(i-c(1))/R);
           f1 = cos(th)*cos(ph) / sqrt(cos(th)^2 + sin(th)^2);
           f2 = cos(th)*sin(ph) / (cos(th)^2 + sin(th)^2);
           f3 = -sin(th) / sqrt(cos(th)^2 + sin(th)^2);       
       else 
           f1 = 0;
           f2 = 0;
           f3 = 0;
       end
%BIPOLAR DROPLET in x direction
     elseif strcmp(dirfield, 'BIPOLARX')
       if r(i,j,k,c) < R 
           th = (i-c(3))/(R*r(i/R,j/R,k/R,c/R));
           ph = atan2((j-c(2))/R,(k-c(3))/R);
           f3 = -cos(th)*cos(ph) / sqrt(cos(th)^2 + sin(th)^2);
           f2 = cos(th)*sin(ph) / (cos(th)^2 + sin(th)^2);
           f1 = -sin(th) / sqrt(cos(th)^2 + sin(th)^2);       
       else 
           f1 = 0;
           f2 = 0;
           f3 = 0;
       end
%read DIR from file   
     elseif strcmp(dirfield, 'FILE')  
         f1 = dirarray(3*ind(i,j,k)-2);
         f2 = dirarray(3*ind(i,j,k)-1);
         f3 = dirarray(3*ind(i,j,k)-0);
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
function [eps,mu] = eblock(i,j,k)
    
    if abs(gain) > 0 && r(i,j,k,c) < R
        G = - 1i * gain * NO;
    else
        G = 0;
    end
    
    %PML factors
    if DPMLs(1) == 0
        sxs = 1;
    elseif i <= abs(DPMLs(1))
        sxs = s(i,nx,-DPMLs(1));
    else
        sxs = 1;
    end
             
    if DPMLs(2) == 0
        sys = 1;
    elseif j <= abs(DPMLs(2))
        sys = s(j,ny,-DPMLs(2));
    else
        sys = 1;
    end
            
    if DPMLs(3) == 0
        szs = 1;
    elseif k <= abs(DPMLs(3))
        szs = s(k,nz,-DPMLs(3));
    else
        szs = 1;    
    end

    if DPMLe(1) == 0
        sxe = 1;
    elseif i > nx - DPMLe(1)
        sxe = s(i,nx,DPMLe(1));
    else
        sxe = 1;    
    end
             
    if DPMLe(2) == 0
        sye = 1;
    elseif j > ny - DPMLe(2)
        sye = s(j,ny,DPMLe(2));
    else
        sye = 1;    
    end
            
    if DPMLe(3) == 0
        sze = 1;
    elseif k > nz - DPMLe(3)
        sze = s(k,nz,DPMLe(3));
    else
        sze = 1;    
    end
    
    
    [dir1,dir2,dir3]=dir(i,j,k);
    Qxx = (dir1*dir1 * 3 -1)/2;
    Qxy = (dir1*dir2 * 3)/2;
    Qxz = (dir1*dir3 * 3)/2;
    Qyy = (dir2*dir2 * 3 -1)/2;
    Qyz = (dir2*dir3 * 3)/2;
    Qzz = (dir3*dir3 * 3 -1)/2;
    
    Qten = [Qxx Qxy Qxz; Qxy Qyy Qyz; Qxz Qyz Qzz];
    
    eps = (deps(i,j,k) * Qten + (epsO(i,j,k)+deps(i,j,k)/3) * eye (3)); 
    eps(1,1) = eps(1,1) * (sys*szs/sxs) * (sye*sze/sxe);
    eps(1,2) = eps(1,2) * (sxs*szs/sys) * (sxe*sze/sye);
    eps(1,3) = eps(1,3) * (sxs*sys/szs) * (sxe*sye/sze);
    eps(2,1) = eps(2,1) * (sys*szs/sxs) * (sye*sze/sxe);
    eps(2,2) = eps(2,2) * (sxs*szs/sys) * (sxe*sze/sye);
    eps(2,3) = eps(2,3) * (sxs*sys/szs) * (sxe*sye/sze);
    eps(3,1) = eps(3,1) * (sys*szs/sxs) * (sye*sze/sxe);
    eps(3,2) = eps(3,2) * (sxs*szs/sys) * (sxe*sze/sye);
    eps(3,3) = eps(3,3) * (sxs*sys/szs) * (sxe*sye/sze);
    

    %mu tensor
    mxx = sqrt((sxs/(sys*szs)) * (sxe/(sye*sze))  );
    myy = sqrt((sys/(sxs*szs)) * (sye/(sxe*sze))  );
    mzz = sqrt((szs/(sxs*sys)) * (sze/(sxe*sye))  );  
    
    %only mu^(-1/2) is needed in calculations 
    mu = [mxx 0 0; 0 myy 0; 0 0 mzz];
    
end

%PML FUNCTIONS
function sout = s(l,n,PML)
    sigma1 = 10*2*lambda2/(4*pi*abs(PML));   
    sout =  psi(l,n,PML) - 1i * sigma1 * eta(l,n,PML);
end

%coordinate stretching
function psiout =  psi(l,n,PML)
     psiout = 1;
%     if PML > 0
%         if l > R
%             psiout = ((n-l)/(n-R))^2;
%         else
%             psiout = 1;
%         end
%     elseif PML < 0
%         if l < R
%             psiout = ((l)/(R))^2;
%         else
%             psiout = 1;
%         end        
%     else
%         psiout = 1;
%     end
end

function etaout = eta(l,n,PML)
    if PML < 0 && l <= abs(PML)
        etaout = (abs((l-abs(PML)))/abs(PML))^3;
    elseif PML > 0 && l > n - PML
        etaout = (abs((l-(n-PML)))/PML)^3;
    else
        etaout = 0;
    end
end
  
function [P2,D2]=sortem(P,D)
    % this function takes in two matrices P and D, presumably the output 
    % from Matlab's eig function, and then sorts the columns of P to 
    % match the sorted columns of D (going from largest to smallest)
    % 
    % EXAMPLE: 
    % 
    % D =
    %    -90     0     0
    %      0   -30     0
    %      0     0   -60
    % P =
    %      1     2     3
    %      1     2     3
    %      1     2     3
    % 
    % [P,D]=sortem(P,D)
    % P =
    %      2     3     1
    %      2     3     1
    %      2     3     1
    % D =
    %    -30     0     0
    %      0   -60     0
    %      0     0   -90
    D2=diag(sort(diag(D),'descend')); % make diagonal matrix out of sorted diagonal values of input D
    [~, ind]=sort(diag(D),'descend'); % store the indices of which columns the sorted eigenvalues come from
    P2=P(:,ind); % arrange the columns in this order
end
        
function A = read_dir(filename)
    fileID = fopen(filename,'r');
    A = fread(fileID,'float');
    fclose(fileID);
end

function [] = exportdir(name)
    dirx = zeros(nx,ny,nz);
    diry = zeros(nx,ny,nz);
    dirz = zeros(nx,ny,nz);
    for di = 1:nx 
        for dj = 1:ny
            for dk = 1:nz
                [dirx(di,dj,dk),diry(di,dj,dk),dirz(di,dj,dk)] = dir(di,dj,dk);
            end
        end
    end
    dirout = cat(4,dirx,diry,dirz);
    dirout = permute(dirout, [4,1,2,3]);
    
    outvec = strcat(num2str(file),'_DIR_',name,'_',num2str(nx),'_',num2str(ny),'_',num2str(nz),'.raw');
    
    fid = fopen(outvec,'w');
    fwrite(fid,dirout,'float');
    fclose(fid);  

end


end