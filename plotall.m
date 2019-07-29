function [] = plotall(Exin,Eyin,Ezin,value,oct,outgen,cycle,Q)
%plots the entire region taking symmetry BC into account

global nx ny nz file c

NX = 2*nx;
NY = 2*ny;
NZ = 2*(nz-1);

Exall = zeros(NX,NY,NZ);
Eyall = zeros(NX,NY,NZ);
Ezall = zeros(NX,NY,NZ);
INT = zeros(NX,NY,NZ);
S0 = zeros(NX,NY,NZ);
S1 = zeros(NX,NY,NZ);
S2 = zeros(NX,NY,NZ);
S3 = zeros(NX,NY,NZ);
S12phi = zeros(NX,NY,NZ);
S12 = zeros(NX,NY,NZ);

if oct == 1
    for i = 1:nx
        for j = 1:ny
            for k = 2:nz
                %1st octant
                Exall(i+nx,j+ny,k+nz) = Exin(i,j,k-1); 
                Eyall(i+nx,j+ny,k+nz) = Eyin(i,j,k-1); 
                Ezall(i+nx,j+ny,k+nz) = Ezin(i,j,k-1); 
                %2nd octant
                Exall(i,j+ny,k+nz) = -Exin(nx-i+1,j,k-1); 
                Eyall(i,j+ny,k+nz) = Eyin(nx-i+1,j,k-1); 
                Ezall(i,j+ny,k+nz) = Ezin(nx-i+1,j,k-1);     
                %3rd octant
                Exall(i,j,k+nz) = -Exin(nx-i+1,ny-j+1,k-1); 
                Eyall(i,j,k+nz) = -Eyin(nx-i+1,ny-j+1,k-1); 
                Ezall(i,j,k+nz) = Ezin(nx-i+1,ny-j+1,k-1); 
                %4th octant
                Exall(i+nx,j,k+nz) = Exin(i,ny-j+1,k-1); 
                Eyall(i+nx,j,k+nz) = -Eyin(i,ny-j+1,k-1); 
                Ezall(i+nx,j,k+nz) = Ezin(i,ny-j+1,k-1);          
                %5th octant
                Exall(i+nx,j+ny,k) = Exin(i,j,nz-k+1-1); 
                Eyall(i+nx,j+ny,k) = Eyin(i,j,nz-k+1-1); 
                Ezall(i+nx,j+ny,k) = -Ezin(i,j,nz-k+1-1); 
                %6th octant
                Exall(i,j+ny,k) = -Exin(nx-i+1,j,nz-k+1-1); 
                Eyall(i,j+ny,k) = Eyin(nx-i+1,j,nz-k+1-1); 
                Ezall(i,j+ny,k) = -Ezin(nx-i+1,j,nz-k+1-1);  
                %7th octant
                Exall(i,j,k) = -Exin(nx-i+1,ny-j+1,nz-k+1-1); 
                Eyall(i,j,k) = -Eyin(nx-i+1,ny-j+1,nz-k+1-1); 
                Ezall(i,j,k) = -Ezin(nx-i+1,ny-j+1,nz-k+1-1); 
                %8th octant
                Exall(i+nx,j,k) = Exin(i,ny-j+1,nz-k+1-1); 
                Eyall(i+nx,j,k) = -Eyin(i,ny-j+1,nz-k+1-1); 
                Ezall(i+nx,j,k) = -Ezin(i,ny-j+1,nz-k+1-1); 
            end
        end
    end
elseif oct == 2
    for i = 1:nx
        for j = 1:ny
            for k = 2:nz
                %1st octant
                Exall(i+nx,j+ny,k+(nz-1)-1) = Exin(nx-i+1,j,k); 
                Eyall(i+nx,j+ny,k+(nz-1)-1) = -Eyin(nx-i+1,j,k); 
                Ezall(i+nx,j+ny,k+(nz-1)-1) = -Ezin(nx-i+1,j,k); 
                %2nd octant
                Exall(i,j,k+(nz-1)-1) = Exin(i,j,k); 
                Eyall(i,j,k+(nz-1)-1) = Eyin(i,j,k); 
                Ezall(i,j,k+(nz-1)-1) = Ezin(i,j,k); 
                %3rd octant
                Exall(i,j,k+(nz-1)-1) = Exin(i,ny-j+1,k); 
                Eyall(i,j,k+(nz-1)-1) = -Eyin(i,ny-j+1,k); 
                Ezall(i,j,k+(nz-1)-1) = Ezin(i,ny-j+1,k); 
                %4th octant
                Exall(i+nx,j,k+(nz-1)-1) = Exin(nx-i+1,ny-j+1,k); 
                Eyall(i+nx,j,k+(nz-1)-1) = Eyin(nx-i+1,ny-j+1,k); 
                Ezall(i+nx,j,k+(nz-1)-1) = -Ezin(nx-i+1,ny-j+1,k); 
                %5th octant
                Exall(i+nx,j+ny,k-1) = Exin(nx-i+1,j,nz-k+1); 
                Eyall(i+nx,j+ny,k-1) = -Eyin(nx-i+1,j,nz-k+1); 
                Ezall(i+nx,j+ny,k-1) = Ezin(nx-i+1,j,nz-k+1); 
                %6th octant
                Exall(i,j+ny,k-1) = Exin(i,j,nz-k+1); 
                Eyall(i,j+ny,k-1) = Eyin(i,j,nz-k+1); 
                Ezall(i,j+ny,k-1) = -Ezin(i,j,nz-k+1);   
                %7th octant
                Exall(i,j,k-1) = Exin(i,ny-j+1,nz-k+1); 
                Eyall(i,j,k-1) = -Eyin(i,ny-j+1,nz-k+1); 
                Ezall(i,j,k-1) = -Ezin(i,ny-j+1,nz-k+1); 
                %8th octant
                Exall(i+nx,j,k-1) = Exin(nx-i+1,ny-j+1,nz-k+1); 
                Eyall(i+nx,j,k-1) = Eyin(nx-i+1,ny-j+1,nz-k+1); 
                Ezall(i+nx,j,k-1) = Ezin(nx-i+1,ny-j+1,nz-k+1); 
            end
        end
    end
elseif oct == 5
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                %1st octant
                Exall(i+nx,j+ny,k+nz) = -Exin(i,j,nz-k+1); 
                Eyall(i+nx,j+ny,k+nz) = -Eyin(i,j,nz-k+1); 
                Ezall(i+nx,j+ny,k+nz) = Ezin(i,j,nz-k+1); 
                %2nd octant
                Exall(i,j+ny,k+nz) = Exin(nx-i+1,j,nz-k+1); 
                Eyall(i,j+ny,k+nz) = -Eyin(nx-i+1,j,nz-k+1); 
                Ezall(i,j+ny,k+nz) = Ezin(nx-i+1,j,nz-k+1);     
                %3rd octant
                Exall(i,j,k+nz) = Exin(nx-i+1,ny-j+1,nz-k+1); 
                Eyall(i,j,k+nz) = Eyin(nx-i+1,ny-j+1,nz-k+1); 
                Ezall(i,j,k+nz) = Ezin(nx-i+1,ny-j+1,nz-k+1); 
                %4th octant
                Exall(i+nx,j,k+nz) = -Exin(i,ny-j+1,nz-k+1); 
                Eyall(i+nx,j,k+nz) = Eyin(i,ny-j+1,nz-k+1); 
                Ezall(i+nx,j,k+nz) = Ezin(i,ny-j+1,nz-k+1);          
                %5th octant
                Exall(i+nx,j+ny,k) = Exin(i,j,k); 
                Eyall(i+nx,j+ny,k) = Eyin(i,j,k); 
                Ezall(i+nx,j+ny,k) = Ezin(i,j,k); 
                %6th octant
                Exall(i,j+ny,k) = -Exin(nx-i+1,j,k); 
                Eyall(i,j+ny,k) = Eyin(nx-i+1,j,k); 
                Ezall(i,j+ny,k) = Ezin(nx-i+1,j,k);  
                %7th octant
                Exall(i,j,k) = -Exin(nx-i+1,ny-j+1,k); 
                Eyall(i,j,k) = -Eyin(nx-i+1,ny-j+1,k); 
                Ezall(i,j,k) = Ezin(nx-i+1,ny-j+1,k); 
                %8th octant
                Exall(i+nx,j,k) = Exin(i,ny-j+1,k); 
                Eyall(i+nx,j,k) = -Eyin(i,ny-j+1,k); 
                Ezall(i+nx,j,k) = Ezin(i,ny-j+1,k); 
            end
        end
    end
    elseif oct == 6
    for i = 1:nx
        for j = 1:ny
            for k = 2:nz
                %1st octant
                Exall(i+nx,j+ny,k+nz-2) = -Exin(nx-i+1,j,nz-k+1); 
                Eyall(i+nx,j+ny,k+nz-2) = Eyin(nx-i+1,j,nz-k+1); 
                Ezall(i+nx,j+ny,k+nz-2) = -Ezin(nx-i+1,j,nz-k+1); 
                %2nd octant
                Exall(i,j+ny,k+nz-2) = -Exin(i,j,nz-k+1); 
                Eyall(i,j+ny,k+nz-2) = -Eyin(i,j,nz-k+1); 
                Ezall(i,j+ny,k+nz-2) = Ezin(i,j,nz-k+1); 
                %3rd octant
                Exall(i,j,k+nz-2) = -Exin(i,ny-j+1,nz-k+1); 
                Eyall(i,j,k+nz-2) = Eyin(i,ny-j+1,nz-k+1); 
                Ezall(i,j,k+nz-2) = Ezin(i,ny-j+1,nz-k+1); 
                %4th octant
                Exall(i+nx,j,k+nz-2) = -Exin(nx-i+1,ny-j+1,nz-k+1); 
                Eyall(i+nx,j,k+nz-2) = -Eyin(nx-i+1,ny-j+1,nz-k+1); 
                Ezall(i+nx,j,k+nz-2) = -Ezin(nx-i+1,ny-j+1,nz-k+1); 
                %5th octant
                Exall(i+nx,j+ny,k-1) = Exin(nx-i+1,j,k); 
                Eyall(i+nx,j+ny,k-1) = -Eyin(nx-i+1,j,k); 
                Ezall(i+nx,j+ny,k-1) = -Ezin(nx-i+1,j,k); 
                %6th octant
                Exall(i,j+ny,k-1) = Exin(i,j,k); 
                Eyall(i,j+ny,k-1) = Eyin(i,j,k); 
                Ezall(i,j+ny,k-1) = Ezin(i,j,k);   
                %7th octant
                Exall(i,j,k-1) = Exin(i,ny-j+1,k); 
                Eyall(i,j,k-1) = -Eyin(i,ny-j+1,k); 
                Ezall(i,j,k-1) = Ezin(i,ny-j+1,k); 
                %8th octant
                Exall(i+nx,j,k-1) = Exin(nx-i+1,ny-j+1,k); 
                Eyall(i+nx,j,k-1) = Eyin(nx-i+1,ny-j+1,k); 
                Ezall(i+nx,j,k-1) = -Ezin(nx-i+1,ny-j+1,k); 
            end
        end
    end
elseif oct == 7
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                %1st octant
                Exall(i+nx,j+ny,k+nz) = Exin(nx-i+1,ny-j+1,nz-k+1); 
                Eyall(i+nx,j+ny,k+nz) = Eyin(nx-i+1,ny-j+1,nz-k+1); 
                Ezall(i+nx,j+ny,k+nz) = Ezin(nx-i+1,ny-j+1,nz-k+1); 
                %2nd octant
                Exall(i,j+ny,k+nz) = Exin(i,ny-j+1,nz-k+1); 
                Eyall(i,j+ny,k+nz) = -Eyin(i,ny-j+1,nz-k+1); 
                Ezall(i,j+ny,k+nz) = -Ezin(i,ny-j+1,nz-k+1); 
                %3rd octant
                Exall(i,j,k+nz) = -Exin(i,j,nz-k+1); 
                Eyall(i,j,k+nz) = -Eyin(i,j,nz-k+1); 
                Ezall(i,j,k+nz) = Ezin(i,j,nz-k+1); 
                %4th octant
                Exall(i+nx,j,k+nz) = -Exin(nx-i+1,j,nz-k+1); 
                Eyall(i+nx,j,k+nz) = Eyin(nx-i+1,j,nz-k+1); 
                Ezall(i+nx,j,k+nz) = -Ezin(nx-i+1,j,nz-k+1); 
                %5th octant
                Exall(i+nx,j+ny,k) = -Exin(nx-i+1,ny-j+1,k); 
                Eyall(i+nx,j+ny,k) = -Eyin(nx-i+1,ny-j+1,k); 
                Ezall(i+nx,j+ny,k) = Ezin(nx-i+1,ny-j+1,k); 
                %6th octant
                Exall(i,j+ny,k) = -Exin(i,ny-j+1,k); 
                Eyall(i,j+ny,k) = Eyin(i,ny-j+1,k); 
                Ezall(i,j+ny,k) = -Ezin(i,ny-j+1,k);   
                %7th octant
                Exall(i,j,k) = Exin(i,j,k); 
                Eyall(i,j,k) = Eyin(i,j,k); 
                Ezall(i,j,k) = Ezin(i,j,k); 
                %8th octant
                Exall(i+nx,j,k) = Exin(nx-i+1,j,k); 
                Eyall(i+nx,j,k) = -Eyin(nx-i+1,j,k); 
                Ezall(i+nx,j,k) = -Ezin(nx-i+1,j,k); 
            end
        end
    end
elseif oct == 8
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                %1st octant
                Exall(i+nx,j+ny,k+nz) = Exin(i,ny-j+1,nz-k+1); 
                Eyall(i+nx,j+ny,k+nz) = -Eyin(i,ny-j+1,nz-k+1); 
                Ezall(i+nx,j+ny,k+nz) = -Ezin(i,ny-j+1,nz-k+1); 
                %2nd octant
                Exall(i,j+ny,k+nz) = -Exin(nx-i+1,ny-j+1,nz-k+1); 
                Eyall(i,j+ny,k+nz) = -Eyin(nx-i+1,ny-j+1,nz-k+1); 
                Ezall(i,j+ny,k+nz) = -Ezin(nx-i+1,ny-j+1,nz-k+1); 
                %3rd octant
                Exall(i,j,k+nz) = Exin(nx-i+1,j,nz-k+1); 
                Eyall(i,j,k+nz) = -Eyin(nx-i+1,j,nz-k+1); 
                Ezall(i,j,k+nz) = Ezin(nx-i+1,j,nz-k+1); 
                %4th octant
                Exall(i+nx,j,k+nz) = -Exin(i,j,nz-k+1); 
                Eyall(i+nx,j,k+nz) = -Eyin(i,j,nz-k+1); 
                Ezall(i+nx,j,k+nz) = Ezin(i,j,nz-k+1); 
                %5th octant
                Exall(i+nx,j+ny,k) = -Exin(i,ny-j+1,k); 
                Eyall(i+nx,j+ny,k) = Eyin(i,ny-j+1,k); 
                Ezall(i+nx,j+ny,k) = -Ezin(i,ny-j+1,k); 
                %6th octant
                Exall(i,j+ny,k) = Exin(nx-i+1,ny-j+1,k); 
                Eyall(i,j+ny,k) = Eyin(nx-i+1,ny-j+1,k); 
                Ezall(i,j+ny,k) = -Ezin(nx-i+1,ny-j+1,k);   
                %7th octant
                Exall(i,j,k) = -Exin(nx-i+1,j,k); 
                Eyall(i,j,k) = Eyin(nx-i+1,j,k); 
                Ezall(i,j,k) = Ezin(nx-i+1,j,k); 
                %8th octant
                Exall(i+nx,j,k) = Exin(i,j,k); 
                Eyall(i+nx,j,k) = Eyin(i,j,k); 
                Ezall(i+nx,j,k) = Ezin(i,j,k); 
            end
        end
    end    
else
    disp('Wrong octant');
end

reEx = real(Exall);
reEy = real(Eyall);
reEz = real(Ezall);

for i = 1:NX
    for j = 1:NY
        for k = 1:NZ
            INT(i,j,k) = reEx(i,j,k)^2 + reEy(i,j,k)^2 + reEz(i,j,k)^2;
            %STOKES PARAMETERS
            S0(i,j,k) = abs(Exall(i,j,k))^2 + abs(Eyall(i,j,k))^2;
            S1(i,j,k) = S0(i,j,k)^(-1) * (abs(Exall(i,j,k))^2 - abs(Eyall(i,j,k))^2);
            S2(i,j,k) = 2*S0(i,j,k)^(-1) * real(conj(Exall(i,j,k))*Eyall(i,j,k));
%             S3(i,j,k) = 2*SO^(-1) * imag(conj(Exall(i,j,k))*Eyall(i,j,k));
            %STOKES FIELD PHASE
            S12phi(i,j,k) = atan2(S2(i,j,k),S1(i,j,k));
            %STOKES FIELD
%             S12(i,j,k) = sqrt(S1(i,j,k)^2 + S2(i,j,k)^2) * exp(1i * S12phi(i,j,k));
            
            
        end
    end
end


Ex = Exall;
Ey = Eyall;
Ez = Ezall;
intE = INT;

clear Exall Eyall Ezall INT;
    
px = fix(nx);
py = fix(ny); 
pz = fix(nz);

%FUNCTIONS TO PLOT E
    sE1 = real(squeeze(Ex(:,:,pz)));
    sE2 = real(squeeze(Ey(:,:,pz)));
    sE3 = real(squeeze(Ez(:,:,pz)));
    sE4 = real(squeeze(Ex(:,py,:)));
    sE5 = real(squeeze(Ey(:,py,:)));
    sE6 = real(squeeze(Ez(:,py,:)));
    sE7 = real(squeeze(Ex(px,:,:)));
    sE8 = real(squeeze(Ey(px,:,:)));
    sE9 = real(squeeze(Ez(px,:,:)));
    lE1 = real(squeeze(Ex(:,py,pz)));
    lE2 = real(squeeze(Ey(:,py,pz)));
    lE3 = real(squeeze(Ez(:,py,pz)));
    lE4 = real(squeeze(Ex(px,:,pz)));
    lE5 = real(squeeze(Ey(px,:,pz)));
    lE6 = real(squeeze(Ez(px,:,pz)));
    lE7 = real(squeeze(Ex(px,py,:)));
    lE8 = real(squeeze(Ey(px,py,:)));
    lE9 = real(squeeze(Ez(px,py,:)));
    
    if c(3) > nz/2
        d = -1;
    else
        d = 1;
    end
    
    S12phi1 = real(squeeze(S12phi(:,:,pz - fix(nz/2))));  
    S12phi2 = real(squeeze(S12phi(:,:,pz)));    
    S12phi3 = real(squeeze(S12phi(:,:,pz + fix(nz/2))));
    
    intE1 = real(squeeze(intE(:,:,pz - fix(nz/2))));  
    intE2 = real(squeeze(intE(:,:,pz)));
    intE3 = real(squeeze(intE(:,:,pz + fix(nz/2))));
    
    intE31 = real(squeeze(intE(:,:,2*nz - fix(pi/(2*sqrt(real(value))))))); 
    sE11 = real(squeeze(Ex(:,:,2*nz - fix(pi/(2*sqrt(real(value)))))));
    sE21 = real(squeeze(Ey(:,:,2*nz - fix(pi/(2*sqrt(real(value)))))));
%     intE31 = real(squeeze(intE(:,:,1))); 
%     sE11 = real(squeeze(Ex(:,:,1)));
%     sE21 = real(squeeze(Ey(:,:,1)));
%     intE1 = real(squeeze(intE(:,:,20)));  
%     intE2 = real(squeeze(intE(:,:,40)));
%     intE3 = real(squeeze(intE(:,:,51)));   

    intE4 = real(squeeze(intE(:,py - fix(ny/2),:)));  
    intE5 = real(squeeze(intE(:,py,:)));
    intE6 = real(squeeze(intE(:,py + fix(ny/2),:)));  
    intE7 = real(squeeze(intE(px - fix(nx/2),:,:)));  
    intE8 = real(squeeze(intE(px,:,:)));
    intE9 = real(squeeze(intE(px + fix(nx/2),:,:)));  
    
    titles1 = ["Ex(x,y)","Ey(x,y)","Ez(x,y)","Ex(x,z)","Ey(x,z)","Ez(x,z)","Ex(y,z)","Ey(y,z)","Ez(y,z)"];
    titles2 = ["Ex(x)","Ey(x)","Ez(x)","Ex(y)","Ey(y)","Ez(y)","Ex(z)","Ey(z)","Ez(z)"];
    titlesintE = ["Eint(z/4)","Eint(z/2)","Eint(3*z/4)","Eint(y/4)","Eint(y/2)","Eint(3y/4)","Eint(x/4)","Eint(x/2)","Eint(3x/4)"];

    clear Ex Ey Ez

    %določi min in max vseh 9 grafov
    maxar = zeros (1,9);
    minar = zeros (1,9);
    maxint = zeros (1,9);
    minint = zeros (1,9);
    for i = 1:9
       maxar(1,i) = max(eval(['sE', num2str(i),'(:)']));
       minar(1,i) = min(eval(['sE', num2str(i),'(:)']));

       maxint(1,i) = max(intE(:));
       minint(1,i) = min(intE(:));   

    end

    %nariše fig1 in subfigure
    clear zlim
    fig5=figure;
    set(fig5,'visible','off');
    for i = 1:9
       subplot(3,3,i);
       surf(eval(['sE', num2str(i)]),'edgecolor','none');
       if 1.1*max(maxar)==0 && 1.1*min(minar)==0
           zlim([-1 1]);
       else
           zlim([1.1*min(minar) 1.1*max(maxar)]);
       end
       caxis([1.1*min(minar) 1.1*max(maxar)]);
       title(titles1(i))
       daspect([1 1 1])
       %shading interp;
       view(2)
    end

    %nastavi pozicijo colormap legende
    hp4 = get(subplot(3,3,6),'Position');
    colorbar('Position', [1.33*hp4(1)  hp4(2)/2  hp4(3)/10  3*hp4(4)])

    %nariše intenziteta
    clear zlim
    fig6=figure;
    set(fig6,'visible','off');
    for i = 1:9
       subplot(3,3,i);
       surf(eval(['intE', num2str(i)]),'edgecolor','none');
       zlim([1.1*min(minint) 1.1*max(maxint)]);
       caxis([1.1*min(minint) 1.1*max(maxint)]);
       title(titlesintE(i))
       daspect([1 1 1])
       view(2)
    end
    
    fig7=figure;
    set(fig7,'visible','off');
    for i = 1:9
       subplot(3,3,i);
       plot(eval(['lE', num2str(i)]));
       if 1.1*min(minar)==0 && 1.1*max(maxar)==0 
           ylim([-1 1]);
       else
           ylim([1.1*min(minar) 1.1*max(maxar)]);
       end
       title(titles2(i))
    end
    
    if abs(Q) > 0
        fig8=figure;
        set(fig8,'visible','off');
        surf(intE2','edgecolor','none');
        zlim([-1.1*max(maxint) 1.1*max(maxint)]);
        caxis([1.1*min(minint) 1.1*max(maxint)]);
        daspect([1 1 1])
        hold on
        [x,y]=meshgrid(1:2*nx,1:2*ny);
        hmany = 3;
        quiver(x(1:hmany:2*ny,1:hmany:2*nx),y(1:hmany:2*ny,1:hmany:2*nx),sE1(1:hmany:2*nx,1:hmany:2*ny)',sE2(1:hmany:2*nx,1:hmany:2*ny)','w');
        hold off
        view(2)
        
        fig9=figure;
        set(fig9,'visible','off');
        surf(intE31','edgecolor','none');
        zlim([-1.1*max(maxint) 1.1*max(maxint)]);
        caxis([1.1*min(minint) 1.1*max(maxint)]);
        daspect([1 1 1])
        hold on
        [x,y]=meshgrid(1:2*nx,1:2*ny);
        hmany = 3;
        quiver(x(1:hmany:2*ny,1:hmany:2*nx),y(1:hmany:2*ny,1:hmany:2*nx),sE11(1:hmany:2*nx,1:hmany:2*ny)',sE21(1:hmany:2*nx,1:hmany:2*ny)','w');
        hold off
        view(2)
    end
    
    fig10=figure;
    set(fig10,'visible','off');
    for i = 1:3
       subplot(1,3,i);
       surf(eval(['S12phi', num2str(i)]),'edgecolor','none');
       zlim([-pi pi]);
       caxis([-pi pi]);
       daspect([1 1 1])
       view(2)
    end
    
    clear sE1 sE2 sE3 sE4 sE5 sE6 sE7 sE8 sE9
    clear lE1 lE2 lE3 lE4 lE5 lE6 lE7 lE8 lE9

    %WHERE TO SAVE THEM
    if isempty(file) == 0
        out = strcat(file,'/plotall/',num2str(real(2*pi/sqrt(value))),'_',num2str(Q),'_',outgen,'_');
        print(fig5,strcat(out,'E2Dplane_n',num2str(cycle),'.png'),'-dpng');
        pause(1) 
        print(fig6,strcat(out,'Eintensty',num2str(cycle),'.png'),'-dpng');
        pause(1) 
        print(fig7,strcat(out,'E1D_cycle',num2str(cycle),'.png'),'-dpng');
        pause(1) 
        if abs(Q) > 0
            print(fig8,strcat(out,'polarization_mid',num2str(cycle),'.png'),'-dpng');
            pause(1) 
            print(fig9,strcat(out,'polarization_end',num2str(cycle),'.png'),'-dpng');
            pause(1) 
        end
        print(fig10,strcat(out,'Stokes_phase',num2str(cycle),'.png'),'-dpng');
    end

clear fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8;


end
