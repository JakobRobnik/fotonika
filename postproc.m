global nx ny nz NX NY NZ file Efile lambda out

file = strcat('/data/Matlab/UPML/analizemodes/');
Efile = strcat(file,'Efield_lda9.2463-0.0042056iQ_1099.2855_40x40x35.mat');
load(Efile);

lambda = 9.2463;
omega = (2*pi)^2/(lambda)^2;

out = strcat(file,'mode_',num2str(lambda),'/');
mkdir(out);

plot = 'ALL';
% plot = 'NOPML';

%box size in px
nx = 40; 
ny = 40;
nz = 35;
NX = nx;
NY = ny;
NZ = nz;
R = 30;
%which octant is simulated
oct = 5;
%thickness of PML, s = at start of domain (i=0) e = end of domain (i=n)
%either or both can be applied
DPMLs = [0,0,0];
DPMLe = [7,7,0];
%center of the droplet
c = [0,0,nz+1];

%calculate components
[reEx,reEy,reEz,imEx,imEy,imEz,intE] = components(vecE,plot);

%if we only simulate one octant - simmetries
[reEx,reEy,reEz] = mirror(reEx,reEy,reEz,oct);
[imEx,imEy,imEz] = mirror(imEx,imEy,imEz,oct);
NX = 2*nx;
NY = 2*ny;
NZ = 2*nz;

%calculates phase between components - NOT USEFUL
%[phx,phy,phz] = phase(reEx,reEy,reEz,imEx,imEy,imEz);

%calculates stokes' fields
[S0,S1,S2,S3,S12ph,S12] = stokes(reEx+1i*imEx,reEy+1i*imEy);

%output for paraview
paraview_out_field(reEx,reEy,reEz,'realE');
paraview_out_field(imEx,imEy,imEz,'imagE');

paraview_out_scalar(phx,'phaseX');
paraview_out_scalar(phy,'phaseY');
paraview_out_scalar(phz,'phaseZ');

paraview_out_scalar(S0,'stokes0');
paraview_out_scalar(S1,'stokes1');
paraview_out_scalar(S2,'stokes2');
paraview_out_scalar(S3,'stokes3');
paraview_out_scalar(S12ph,'stokes12ph');
paraview_out_scalar(S12,'stokes12');

[reExrot,reEyrot,reEzrot] = rotate(reEx,reEy,reEz,pi/4);
%tukaj dodaj sestevanje 2 vektorjev
Extot = reExrot - reEx;
Eytot = reEyrot - reEy;
Eztot = reEzrot - reEz;

% Extot = zeros(NX,NY,NZ);
% Eytot = zeros(NX,NY,NZ);
% Eztot = zeros(NX,NY,NZ);

% for i = 1:NX
%     for j = 1:NY
%         for k = 1:NZ
%             Extot(i,j,k) = reExrot(i,j,k)+reEx(i,j,k);
%             Eytot(i,j,k) = reEyrot(i,j,k)+reEy(i,j,k);        
%             Eztot(i,j,k) = reEzrot(i,j,k)+reEz(i,j,k);
%         end
%     end
% end

intEtot = zeros(NX,NY,NZ);
intErot = zeros(NX,NY,NZ);

for i = 1:NX
    for j = 1:NY
        for k = 1:NZ
            intE(i,j,k) = (reEx(i,j,k))^2+(reEy(i,j,k))^2+(reEz(i,j,k))^2;
            intErot(i,j,k) = (reExrot(i,j,k))^2+(reEyrot(i,j,k))^2+(reEzrot(i,j,k))^2;            
            intEtot(i,j,k) = (Extot(i,j,k))^2+(Eytot(i,j,k))^2+(Eztot(i,j,k))^2;
        end
    end
end

plotfield(reEx,reEy,reEz,intE,NX/2,NY/2,NZ/2,'original')
plotfield(reExrot,reEyrot,reEzrot,intErot,NX/2,NY/2,NZ/2,'rotated')
plotfield(Extot,Eytot,Eztot,intEtot,NX/2,NY/2,NZ/2,'add_rotated')




function [reEx,reEy,reEz,imEx,imEy,imEz,intE] = components(vecE,plot)
    global nx ny nz
    reEx = zeros(nx,ny,nz);
    reEy = zeros(nx,ny,nz);
    reEz = zeros(nx,ny,nz);
    imEx = zeros(nx,ny,nz);
    imEy = zeros(nx,ny,nz);
    imEz = zeros(nx,ny,nz);
    intE = zeros(nx,ny,nz);
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if strcmp(plot, 'ALL')      
                    reEx(i,j,k) = real(vecE((ind(i,j,k)),1));
                    reEy(i,j,k) = real(vecE((ind(i,j,k))+ nx*ny*nz,1));
                    reEz(i,j,k) = real(vecE((ind(i,j,k))+ 2*nx*ny*nz,1));
                    imEx(i,j,k) = imag(vecE((ind(i,j,k)),1));
                    imEy(i,j,k) = imag(vecE((ind(i,j,k))+ nx*ny*nz,1));
                    imEz(i,j,k) = imag(vecE((ind(i,j,k))+ 2*nx*ny*nz,1));
                    intE(i,j,k) = reEx(i,j,k)^2 + reEy(i,j,k)^2+ reEz(i,j,k)^2;
                elseif strcmp(plot, 'NOPML')
                    if i > nx - DPMLe(1) || j > ny - DPMLe(2) || k > nz - DPMLe(3) || i < DPMLs(1) || j < DPMLs(2) || k < DPMLs(3)
                        reEx(i,j,k) = 0;
                        reEy(i,j,k) = 0;
                        reEz(i,j,k) = 0;
                        imEx(i,j,k) = 0;
                        imEy(i,j,k) = 0;
                        imEz(i,j,k) = 0;  
                    else
                        reEx(i,j,k) = real(vecE((ind(i,j,k)),1));
                        reEy(i,j,k) = real(vecE((ind(i,j,k))+ nx*ny*nz,1));
                        reEz(i,j,k) = real(vecE((ind(i,j,k))+ 2*nx*ny*nz,1));
                        imEx(i,j,k) = imag(vecE((ind(i,j,k)),1));
                        imEy(i,j,k) = imag(vecE((ind(i,j,k))+ nx*ny*nz,1));
                        imEz(i,j,k) = imag(vecE((ind(i,j,k))+ 2*nx*ny*nz,1));                    
                    end
                    intE(i,j,k) = reEx(i,j,k)^2 + reEy(i,j,k)^2 + reEz(i,j,k)^2;
                else
                    disp('Value of plot input argument must be ALL or NOPML');
                end
            end
        end
    end
end

function f = ind(i,j,k)
    global NX NY
    f = (k-1)*NX*NY + (j-1)*NX + (i-1) + 1;
end

function [Exall,Eyall,Ezall] = mirror(Exin,Eyin,Ezin,oct)
    global nx ny nz
    NX = 2*nx;
    NY = 2*ny;
    NZ = 2*nz;

    Exall = zeros(NX,NY,NZ);
    Eyall = zeros(NX,NY,NZ);
    Ezall = zeros(NX,NY,NZ);

    if oct == 1
        for i = 1:nx
            for j = 1:ny
                for k = 1:nz
                    %1st octant
                    Exall(i+nx,j+ny,k+nz) = Exin(i,j,k); 
                    Eyall(i+nx,j+ny,k+nz) = Eyin(i,j,k); 
                    Ezall(i+nx,j+ny,k+nz) = Ezin(i,j,k); 
                    %2nd octant
                    Exall(i,j+ny,k+nz) = -Exin(nx-i+1,j,k); 
                    Eyall(i,j+ny,k+nz) = Eyin(nx-i+1,j,k); 
                    Ezall(i,j+ny,k+nz) = Ezin(nx-i+1,j,k);     
                    %3rd octant
                    Exall(i,j,k+nz) = -Exin(nx-i+1,ny-j+1,k); 
                    Eyall(i,j,k+nz) = -Eyin(nx-i+1,ny-j+1,k); 
                    Ezall(i,j,k+nz) = Ezin(nx-i+1,ny-j+1,k); 
                    %4th octant
                    Exall(i+nx,j,k+nz) = Exin(i,ny-j+1,k); 
                    Eyall(i+nx,j,k+nz) = -Eyin(i,ny-j+1,k); 
                    Ezall(i+nx,j,k+nz) = Ezin(i,ny-j+1,k);          
                    %5th octant
                    Exall(i+nx,j+ny,k) = Exin(i,j,nz-k+1); 
                    Eyall(i+nx,j+ny,k) = Eyin(i,j,nz-k+1); 
                    Ezall(i+nx,j+ny,k) = -Ezin(i,j,nz-k+1); 
                    %6th octant
                    Exall(i,j+ny,k) = -Exin(nx-i+1,j,nz-k+1); 
                    Eyall(i,j+ny,k) = Eyin(nx-i+1,j,nz-k+1); 
                    Ezall(i,j+ny,k) = -Ezin(nx-i+1,j,nz-k+1);  
                    %7th octant
                    Exall(i,j,k) = -Exin(nx-i+1,ny-j+1,nz-k+1); 
                    Eyall(i,j,k) = -Eyin(nx-i+1,ny-j+1,nz-k+1); 
                    Ezall(i,j,k) = -Ezin(nx-i+1,ny-j+1,nz-k+1); 
                    %8th octant
                    Exall(i+nx,j,k) = Exin(i,ny-j+1,nz-k+1); 
                    Eyall(i+nx,j,k) = -Eyin(i,ny-j+1,nz-k+1); 
                    Ezall(i+nx,j,k) = -Ezin(i,ny-j+1,nz-k+1); 
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
                for k = 1:nz
                    %1st octant
                    Exall(i+nx,j+ny,k+nz) = -Exin(nx-i+1,j,nz-k+1); 
                    Eyall(i+nx,j+ny,k+nz) = Eyin(nx-i+1,j,nz-k+1); 
                    Ezall(i+nx,j+ny,k+nz) = -Ezin(nx-i+1,j,nz-k+1); 
                    %2nd octant
                    Exall(i,j+ny,k+nz) = -Exin(i,j,nz-k+1); 
                    Eyall(i,j+ny,k+nz) = -Eyin(i,j,nz-k+1); 
                    Ezall(i,j+ny,k+nz) = Ezin(i,j,nz-k+1); 
                    %3rd octant
                    Exall(i,j,k+nz) = -Exin(i,ny-j+1,nz-k+1); 
                    Eyall(i,j,k+nz) = Eyin(i,ny-j+1,nz-k+1); 
                    Ezall(i,j,k+nz) = Ezin(i,ny-j+1,nz-k+1); 
                    %4th octant
                    Exall(i+nx,j,k+nz) = -Exin(nx-i+1,ny-j+1,nz-k+1); 
                    Eyall(i+nx,j,k+nz) = -Eyin(nx-i+1,ny-j+1,nz-k+1); 
                    Ezall(i+nx,j,k+nz) = -Ezin(nx-i+1,ny-j+1,nz-k+1); 
                    %5th octant
                    Exall(i+nx,j+ny,k) = Exin(nx-i+1,j,k); 
                    Eyall(i+nx,j+ny,k) = -Eyin(nx-i+1,j,k); 
                    Ezall(i+nx,j+ny,k) = -Ezin(nx-i+1,j,k); 
                    %6th octant
                    Exall(i,j+ny,k) = Exin(i,j,k); 
                    Eyall(i,j+ny,k) = Eyin(i,j,k); 
                    Ezall(i,j+ny,k) = Ezin(i,j,k);   
                    %7th octant
                    Exall(i,j,k) = Exin(i,ny-j+1,k); 
                    Eyall(i,j,k) = -Eyin(i,ny-j+1,k); 
                    Ezall(i,j,k) = Ezin(i,ny-j+1,k); 
                    %8th octant
                    Exall(i+nx,j,k) = Exin(nx-i+1,ny-j+1,k); 
                    Eyall(i+nx,j,k) = Eyin(nx-i+1,ny-j+1,k); 
                    Ezall(i+nx,j,k) = -Ezin(nx-i+1,ny-j+1,k); 
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
end

function [phx,phy,phz] = phase(reEx,reEy,reEz,imEx,imEy,imEz)
    global NX NY NZ
    phx = zeros(NX,NY,NZ);
    phy = zeros(NX,NY,NZ);
    phz = zeros(NX,NY,NZ);
    for i = 1:NX
        for j = 1:NY
            for k = 1:NZ
                phx(i,j,k) = atan2(imEx(i,j,k),reEx(i,j,k));
                phy(i,j,k) = atan2(imEy(i,j,k),reEy(i,j,k));
                phz(i,j,k) = atan2(imEz(i,j,k),reEz(i,j,k));
            end
        end
    end
end

function [] = paraview_out_scalar(S,name)
    global NX NY NZ out lambda

    nx = NX;
    ny = NY;
    nz = NZ;

    outS = strcat(num2str(out),num2str(name),'_',num2str(nx),'_',num2str(ny),'_',num2str(nz),'_lda',num2str(lambda),'.raw');
    
    fid = fopen(outS,'w');
    fwrite(fid,S,'float');
    fclose(fid);    
end

function [S0,S1,S2,S3,S12ph,S12] = stokes(Vx,Vy)
    global NX NY NZ
    S0 = zeros(NX,NY,NZ);
    S1 = zeros(NX,NY,NZ);
    S2 = zeros(NX,NY,NZ);
    S3 = zeros(NX,NY,NZ);
    S12ph = zeros(NX,NY,NZ);
    S12 = zeros(NX,NY,NZ);
    for i = 1:NX
        for j = 1:NY
            for k = 1:NZ
                %STOKES PARAMETERS
                S0(i,j,k) = abs(Vx(i,j,k))^2 + abs(Vy(i,j,k))^2;
                S1(i,j,k) = S0(i,j,k)^(-1) * (abs(Vx(i,j,k))^2 - abs(Vy(i,j,k))^2);
                S2(i,j,k) = 2*S0(i,j,k)^(-1) * real(conj(Vx(i,j,k))*Vy(i,j,k));
                S3(i,j,k) = 2*S0(i,j,k)^(-1) * imag(conj(Vx(i,j,k))*Vy(i,j,k));
                %STOKES FIELD PHASE
                S12ph(i,j,k) = atan2(S2(i,j,k),S1(i,j,k));
                %STOKES FIELD
                S12(i,j,k) = sqrt(S1(i,j,k)^2 + S2(i,j,k)^2) * exp(1i * S12ph(i,j,k));
            end
        end
    end
end

function [Vx2,Vy2,Vz2] = rotate(Vx,Vy,Vz,angle)
    global NX NY NZ
    Vx2 = zeros(NX,NY,NZ);
    Vy2 = zeros(NX,NY,NZ);
    Vz2 = zeros(NX,NY,NZ);
    for i = 1:NX
        for j = 1:NY
            for k = 1:NZ
                if r(i,j,0,[NX/2,NY/2,0]) < min([NX/2 NY/2 NZ/2])
                    phi = atan2((j-NY/2),(i-NX/2));
                    i2 = round(cos(phi-angle)*sqrt((i-NX/2)^2+(j-NY/2)^2)+NX/2);
                    j2 = round(sin(phi-angle)*sqrt((i-NX/2)^2+(j-NY/2)^2)+NY/2);
                    Vx2(i,j,k) = Vx(i2,j2,k)*cos(angle)-Vy(i2,j2,k)*sin(angle);
                    Vy2(i,j,k) = Vx(i2,j2,k)*sin(angle)+Vy(i2,j2,k)*cos(angle);
                    Vz2(i,j,k) = Vz(i2,j2,k);  
                else
                    Vx2(i,j,k) = 0;
                    Vy2(i,j,k) = 0;
                    Vz2(i,j,k) = 0; 
                end
                
            end
        end
    end
end

function f = r(i,j,k,c)
       f = sqrt((i-c(1))^2 + (j-c(2))^2 + (k-c(3))^2);
end

function [] = plotfield(Ex,Ey,Ez,intE,px,py,pz,name)
    global file lambda NZ NY NX out

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
    
    intE2 = real(squeeze(intE(:,:,pz)));
    
    titles1 = ["Ex(x,y)","Ey(x,y)","Ez(x,y)","Ex(x,z)","Ey(x,z)","Ez(x,z)","Ex(y,z)","Ey(y,z)","Ez(y,z)"];
    
    intE31 = real(squeeze(intE(:,:,NZ - fix(lambda/4)))); 
    sE11 = real(squeeze(Ex(:,:,NZ - fix(lambda/4))));
    sE21 = real(squeeze(Ey(:,:,NZ - fix(lambda/4))));
    
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
    
    fig8=figure;
    set(fig8,'visible','off');
    surf(intE2','edgecolor','none');
    zlim([-1.1*max(maxint) 1.1*max(maxint)]);
    caxis([1.1*min(minint) 1.1*max(maxint)]);
    daspect([1 1 1])
    hold on
    [x,y]=meshgrid(1:NX,1:NY);
    hmany = 3;
    quiver(x(1:hmany:NY,1:hmany:NX),y(1:hmany:NY,1:hmany:NX),sE1(1:hmany:NX,1:hmany:NY)',sE2(1:hmany:NX,1:hmany:NY)','w');
    hold off
    view(2)
        
    fig9=figure;
    set(fig9,'visible','off');
    surf(intE31','edgecolor','none');
    zlim([-1.1*max(maxint) 1.1*max(maxint)]);
    caxis([1.1*min(minint) 1.1*max(maxint)]);
    daspect([1 1 1])
    hold on
    [x,y]=meshgrid(1:NX,1:NY);
    hmany = 3;
    quiver(x(1:hmany:NY,1:hmany:NX),y(1:hmany:NY,1:hmany:NX),sE11(1:hmany:NX,1:hmany:NY)',sE21(1:hmany:NX,1:hmany:NY)','w');
    hold off
    view(2)
    
    if isempty(file) == 0
        out_name = strcat(out,name,'/');
        mkdir(out_name);
        print(fig5,strcat(out_name,'E2Dplane_n.png'),'-dpng');
        pause(1) 
        print(fig8,strcat(out_name,'polarization_mid.png'),'-dpng');
        pause(1) 
        print(fig9,strcat(out_name,'polarization_end.png'),'-dpng');
        pause(1) 
    end
    
end