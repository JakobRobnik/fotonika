 function [odivH,odHmax,odivE,odEmax,dW,mH,mE,type] = script_func_block(vecH,vecE,value,whichplot,cycle,outgen)
%SCRIPT_FUNC calculates and draws the fields

global nx ny nz dx dy dz Exx Exy Exz Eyy Eyz Ezz file pfield R DPMLs DPMLe c;

if c(1)==0 || c(1)==0.5
    di = 1;
else
    di = -R/4;
end
if c(2)==0 || c(2)==0.5
    dj = 1;
else
    dj = -R/4;
end
if c(3)==0 || c(3)==0.5
    dk = 1;
else
    dk = -R/4;
end
    
%which plane to draw
px = int16(nx/2);
py = int16(ny/2);
pz = int16(nz/2);

% px = fix(nx/4);
% py = fix(nx/4); 
% pz = fix(nx/4); %plane in z which we want to show


Hx = zeros(nx,ny,nz);
Hy = zeros(nx,ny,nz);
Hz = zeros(nx,ny,nz);
reEx = zeros(nx,ny,nz);
reEy = zeros(nx,ny,nz);
reEz = zeros(nx,ny,nz);
imEx = zeros(nx,ny,nz);
imEy = zeros(nx,ny,nz);
imEz = zeros(nx,ny,nz);
intH = zeros(nx,ny,nz);
intE = zeros(nx,ny,nz);
divH = 0;
divE = 0;
maxdivH = 0;
maxdivE = 0;
W = 0;
Win = 0;

%REARRANGE THE MATRIX TO FIND COMPONENTS OF H, CALCULATE DIVHs 
for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            if strcmp(whichplot, 'INSIDE')
                if r(i,j,k) > R 
                    Hx(i,j,k) = 0;
                    Hy(i,j,k) = 0;
                    Hz(i,j,k) = 0; 
                    intH(i,j,k) = 0; 
                else
                    Hx(i,j,k) = vecH(3*(ind(i,j,k)) - 2,1);
                    Hy(i,j,k) = vecH(3*(ind(i,j,k)) - 1,1);
                    Hz(i,j,k) = vecH(3*(ind(i,j,k)) - 0,1);
                    intH(i,j,k) = Hx(i,j,k)^2 + Hy(i,j,k)^2 + Hz(i,j,k)^2;
                end
            elseif strcmp(whichplot, 'ALL')
                Hx(i,j,k) = vecH((ind(i,j,k)),1);% * psix * psiy * psiz;
                Hy(i,j,k) = vecH((ind(i,j,k))+ nx*ny*nz,1);% * psix * psiy * psiz;
                Hz(i,j,k) = vecH((ind(i,j,k))+ 2*nx*ny*nz,1);% * psix * psiy * psiz;
                intH(i,j,k) = Hx(i,j,k)^2 + Hy(i,j,k)^2 + Hz(i,j,k)^2;
            elseif strcmp(whichplot, 'GUESS')
                Hx(i,j,k) = vecH((ind(i,j,k)),1);
                Hy(i,j,k) = vecH((ind(i,j,k))+ nx*ny*nz,1);
                Hz(i,j,k) = vecH((ind(i,j,k))+ 2*nx*ny*nz,1);
                intH(i,j,k) = Hx(i,j,k)^2 + Hy(i,j,k)^2+ Hz(i,j,k)^2;
            elseif strcmp(whichplot, 'NOPML')
                if i > nx - DPMLe(1) || j > ny - DPMLe(2) || k > nz - DPMLe(3) || i < DPMLs(1) || j < DPMLs(2) || k < DPMLe(3)
                    Hx(i,j,k) = 0;
                    Hy(i,j,k) = 0;
                    Hz(i,j,k) = 0;
                else
                    Hx(i,j,k) = vecH((ind(i,j,k)),1);% * psix * psiy * psiz;
                    Hy(i,j,k) = vecH((ind(i,j,k))+ nx*ny*nz,1);% * psix * psiy * psiz;
                    Hz(i,j,k) = vecH((ind(i,j,k))+ 2*nx*ny*nz,1);% * psix * psiy * psiz;
                end
                intH(i,j,k) = Hx(i,j,k)^2 + Hy(i,j,k)^2 + Hz(i,j,k)^2;
            else
                disp('Value of whichplot input argument is not correct.');
            end
            divH =+ dHx(Hx,i,j,k) + dHy(Hy,i,j,k) + dHz(Hz,i,j,k);
            dHmax = abs((dHx(Hx,i,j,k) + dHy(Hy,i,j,k) + dHz(Hz,i,j,k))/(sqrt(Hx(i,j,k)^2 + Hy(i,j,k)^2 + Hz(i,j,k)^2)));
            if dHmax > maxdivH
                maxdivH = dHmax;
            end          
        end
    end
end

for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            if strcmp(whichplot, 'INSIDE')
                if r(i,j,k) > R 
                    reEx(i,j,k) = 0;
                    reEy(i,j,k) = 0;
                    reEz(i,j,k) = 0; 
                    imEx(i,j,k) = 0;
                    imEy(i,j,k) = 0;
                    imEz(i,j,k) = 0; 
                    intE(i,j,k) = 0; 
                else
                    reEx(i,j,k) = real(vecE(3*(ind(i,j,k)) - 2,1));
                    reEy(i,j,k) = real(vecE(3*(ind(i,j,k)) - 1,1));
                    reEz(i,j,k) = real(vecE(3*(ind(i,j,k)) - 0,1));
                    imEx(i,j,k) = imag(vecE(3*(ind(i,j,k)) - 2,1));
                    imEy(i,j,k) = imag(vecE(3*(ind(i,j,k)) - 1,1));
                    imEz(i,j,k) = imag(vecE(3*(ind(i,j,k)) - 0,1));
                    intE(i,j,k) = reEx(i,j,k)^2 + reEy(i,j,k)^2 + reEz(i,j,k)^2;
                end
            elseif strcmp(whichplot, 'ALL')      
                
                reEx(i,j,k) = real(vecE((ind(i,j,k)),1));
                reEy(i,j,k) = real(vecE((ind(i,j,k))+ nx*ny*nz,1));
                reEz(i,j,k) = real(vecE((ind(i,j,k))+ 2*nx*ny*nz,1));
                imEx(i,j,k) = imag(vecE((ind(i,j,k)),1));
                imEy(i,j,k) = imag(vecE((ind(i,j,k))+ nx*ny*nz,1));
                imEz(i,j,k) = imag(vecE((ind(i,j,k))+ 2*nx*ny*nz,1));
                intE(i,j,k) = reEx(i,j,k)^2 + reEy(i,j,k)^2+ reEz(i,j,k)^2;
                if r_qrt(i,j,k) < R
                    Win =+ Ex(i,j,k)^2 + Ey(i,j,k)^2+ Ez(i,j,k)^2;
                end
                W =+ Ex(i,j,k)^2 + Ey(i,j,k)^2+ Ez(i,j,k)^2;
            elseif strcmp(whichplot, 'NOPML')
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
                if r_qrt(i,j,k) < R
                    Win =+ reEx(i,j,k)^2 + reEy(i,j,k)^2+ reEz(i,j,k)^2;
                end
                W =+ reEx(i,j,k)^2 + reEy(i,j,k)^2+ reEz(i,j,k)^2;  
            elseif strcmp(whichplot, 'GUESS')
                reEx(i,j,k) = real(vecE((ind(i,j,k)),1));
                reEy(i,j,k) = real(vecE((ind(i,j,k))+ nx*ny*nz,1));
                reEz(i,j,k) = real(vecE((ind(i,j,k))+ 2*nx*ny*nz,1));
                imEx(i,j,k) = imag(vecE((ind(i,j,k)),1));
                imEy(i,j,k) = imag(vecE((ind(i,j,k))+ nx*ny*nz,1));
                imEz(i,j,k) = imag(vecE((ind(i,j,k))+ 2*nx*ny*nz,1));               
                intE(i,j,k) = reEx(i,j,k)^2 + reEy(i,j,k)^2+ reEz(i,j,k)^2;
            else
                disp('Value of whichplot input argument must be ALL or INSIDE');
            end
            divE =+ dHx(reEx,i,j,k) +dHy(reEy,i,j,k) + dHz(reEz,i,j,k);
            dEmax = abs((dHx(reEx,i,j,k) + dHy(reEy,i,j,k) + dHz(reEz,i,j,k))/(sqrt(reEx(i,j,k)^2 + reEy(i,j,k)^2 + reEz(i,j,k)^2)));
            if dEmax > maxdivE
                maxdivE = dEmax;
            end          
        end
    end
end


%PLOTS ENTIRE AREA IF ONLY ONE QUADRANT IS SIMULATED
% if abs(Q) > 0
%     plotall(reEx+4i*imEx,reEy+1i*imEy,reEz+1i*imEz,value,oct,outgen,cycle,Q);
% end

dW = Win/W;
odivH = divH;
odHmax = maxdivH;

mH = [max(max(max(abs(Hx)))), max(max(max(abs(Hy)))), max(max(max(abs(Hz))))];
mE = [max(max(max(abs(reEx)))), max(max(max(abs(reEy)))), max(max(max(abs(reEz))))];

if 100*mH(1) < mH(3)
    type = 'TE';
elseif 100*mH(3) < mH(1)
    type = 'TM';
else
    type = 'UNDEFINED';
end
    
    

odivE = divE;
odEmax = maxdivE;

if strcmp(pfield, 'YES')
   % FUNCTIONS TO PLOT H
    s1 = abs(squeeze(Hx(:,:,pz)));
    s2 = abs(squeeze(Hy(:,:,pz)));
    s3 = abs(squeeze(Hz(:,:,pz)));
    s4 = abs(squeeze(Hx(:,py,:)));
    s5 = abs(squeeze(Hy(:,py,:)));
    s6 = abs(squeeze(Hz(:,py,:)));
    s7 = abs(squeeze(Hx(px,:,:)));
    s8 = abs(squeeze(Hy(px,:,:)));
    s9 = abs(squeeze(Hz(px,:,:)));
    l1 = abs(squeeze(Hx(:,py,pz)));
    l2 = abs(squeeze(Hy(:,py,pz)));
    l3 = abs(squeeze(Hz(:,py,pz)));
    l4 = abs(squeeze(Hx(px,:,pz)));
    l5 = abs(squeeze(Hy(px,:,pz)));
    l6 = abs(squeeze(Hz(px,:,pz)));
    l7 = abs(squeeze(Hx(px,py,:)));
    l8 = abs(squeeze(Hy(px,py,:)));
    l9 = abs(squeeze(Hz(px,py,:)));
    
    if c(3) > nz/2
        d = -1;
    else
        d = 1;
    end
    
    intH1 = abs(squeeze(intH(:,:,pz)));  
    intH2 = abs(squeeze(intH(:,:,pz )));
    intH3 = abs(squeeze(intH(:,:,pz )));
    
    
    titles1 = ["Hx(x,y)","Hy(x,y)","Hz(x,y)","Hx(x,z)","Hy(x,z)","Hz(x,z)","Hx(y,z)","Hy(y,z)","Hz(y,z)"];
    titles2 = ["Hx(x)","Hy(x)","Hz(x)","Hx(y)","Hy(y)","Hz(y)","Hx(z)","Hy(z)","Hz(z)"];
    titlesintH = ["intH(z/2)","intH(z/2 + R/2)","intH(z/2 + R)"];

    clear Hx Hy Hz

    %doloci min in max vseh 9 grafov
    maxar = zeros (1,9);
    minar = zeros (1,9);
    maxint = zeros (1,3);
    minint = zeros (1,3);
    for i = 1:9
       maxar(i) = max(max(eval(['s', num2str(i)])));
       minar(i) = min(min(eval(['s', num2str(i)])));
       if i<=3
            maxint(i) = max(max(eval(['intH', num2str(i)])));
            minint(i) = min(min(eval(['intH', num2str(i)])));   
       end
    end

    %narise fig1 in subfigure
    clear zlim
    fig1=figure;
    set(fig1,'visible','off');
    for i = 1:9
       subplot(3,3,i);
       surf(eval(['s', num2str(i)]));
       zlim([1.1*min(minar) 1.1*max(maxar)]);
       caxis([1.1*min(minar) 1.1*max(maxar)]);
       title(titles1(i))
    end

    %nastavi pozicijo colormap legende
    hp4 = get(subplot(3,3,6),'Position');
    colorbar('Position', [1.33*hp4(1)  hp4(2)/2  hp4(3)/10  3*hp4(4)])

    %narise fig1 in subfigure
    clear zlim
    fig2=figure;
    set(fig2,'visible','off');
    for i = 1:9
       subplot(3,3,i);
       surf(eval(['s', num2str(i)]),'edgecolor','none');
       zlim([1.1*min(minar) 1.1*max(maxar)])
       caxis([1.1*min(minar) 1.1*max(maxar)])
       title(titles1(i))
       daspect([1 1 1])
       view(2)
    end

    %nastavi pozicijo colormap legende
    hp4 = get(subplot(3,3,6),'Position');
    colorbar('Position', [1.33*hp4(1)  hp4(2)/2  hp4(3)/10  3*hp4(4)])



    %narise intenziteto
    clear zlim
    fig3=figure;
    set(fig3,'visible','off');
    for i = 1:3
       subplot(1,3,i);
       surf(eval(['intH', num2str(i)]),'edgecolor','none');
       zlim([1.1*min(minint) 1.1*max(maxint)]);
       caxis([1.1*min(minint) 1.1*max(maxint)]);
       title(titlesintH(i))
       daspect([1 1 1])
       %shading interp;
       view(2)
    end

    clear s1 s2 s3 s4 s5 s6 s7 s8 s9
    clear l1 l2 l3 l4 l5 l6 l7 l8 l9

    %FUNCTIONS TO PLOT E
    sE1 = squeeze(reEx(:,:,pz));
    sE2 = squeeze(reEy(:,:,pz));
    sE3 = squeeze(reEz(:,:,pz));
    sE4 = squeeze(reEx(:,py,:));
    sE5 = squeeze(reEy(:,py,:));
    sE6 = squeeze(reEz(:,py,:));
    sE7 = squeeze(reEx(px,:,:));
    sE8 = squeeze(reEy(px,:,:));
    sE9 = squeeze(reEz(px,:,:));
    lE1 = squeeze(reEx(:,py,pz));
    lE2 = squeeze(reEy(:,py,pz));
    lE3 = squeeze(reEz(:,py,pz));
    lE4 = squeeze(reEx(px,:,pz));
    lE5 = squeeze(reEy(px,:,pz));
    lE6 = squeeze(reEz(px,:,pz));
    lE7 = squeeze(reEx(px,py,:));
    lE8 = squeeze(reEy(px,py,:));
    lE9 = squeeze(reEz(px,py,:));
    
    if c(3) > nz/2
        d = -1;
    else
        d = 1;
    end
    
    
    intE1 = abs(squeeze(intE(:,:,pz)));  
    intE2 = abs(squeeze(intE(:,:,pz )));
    intE3 = abs(squeeze(intE(:,:,pz )));   
    intE4 = abs(squeeze(intE(:,pz,:)));  
    intE5 = abs(squeeze(intE(:,pz ,:)));
    intE6 = abs(squeeze(intE(:,pz ,:)));  
    intE7 = abs(squeeze(intE(pz,:,:)));  
    intE8 = abs(squeeze(intE(pz ,:,:)));
    intE9 = abs(squeeze(intE(pz ,:,:)));  
    
    titles1 = ["Ex(x,y)","Ey(x,y)","Ez(x,y)","Ex(x,z)","Ey(x,z)","Ez(x,z)","Ex(y,z)","Ey(y,z)","Ez(y,z)"];
    titles2 = ["Ex(x)","Ey(x)","Ez(x)","Ex(y)","Ey(y)","Ez(y)","Ex(z)","Ey(z)","Ez(z)"];
    titlesintE = ["Eint(z/2)","Eint(z/2 + R/2)","Eint(z/2 + R)","Eint(y/2)","Eint(y/2 + R/2)","Eint(y/2 + R)","Eint(x/2)","Eint(x/2 + R/2)","Eint(x/2 + R)"];

    clear Ex Ey Ez

    %doloci min in max vseh 9 grafov
    maxar = zeros (1,9);
    minar = zeros (1,9);
    maxint = zeros (1,9);
    minint = zeros (1,9);
    for i = 1:9
    %     if mod(i,3) == 0
    %         continue;
    %     end
       maxar(1,i) = max(eval(['sE', num2str(i),'(:)']));
       minar(1,i) = min(eval(['sE', num2str(i),'(:)']));

%        maxint(i) = max(max(eval(['intE', num2str(i)])));
%        minint(i) = min(min(eval(['intE', num2str(i)]))); 
       maxint(1,i) = max(intE(:));
       minint(1,i) = min(intE(:));   

    end

%     rise odvisnost v eni smeri
%     fig3=figure;
%     set(fig3,'visible','off');
%     for i = 1:9
%        subplot(3,3,i);
%        plot(eval(['lE', num2str(i)]));
%        ylim([1.1*min(minar) 1.1*max(maxar)])
%        title(titles2(i))
%     end
    
    %narise fig1 in subfigure
    clear zlim
    fig4=figure;
    set(fig4,'visible','off');
    for i = 1:9
       subplot(3,3,i);
       surf(eval(['sE', num2str(i)]));
       zlim([1.1*min(minar) 1.1*max(maxar)]);
       caxis([1.1*min(minar) 1.1*max(maxar)]);
       title(titles1(i))
    end

    %nastavi pozicijo colormap legende
    hp4 = get(subplot(3,3,6),'Position');
    colorbar('Position', [1.33*hp4(1)  hp4(2)/2  hp4(3)/10  3*hp4(4)])

    %narise fig1 in subfigure
    clear zlim
    fig5=figure;
    set(fig5,'visible','off');
    for i = 1:9
       subplot(3,3,i);
       surf(eval(['sE', num2str(i)]),'edgecolor','none');
       zlim([1.1*min(minar) 1.1*max(maxar)]);
       caxis([1.1*min(minar) 1.1*max(maxar)]);
       title(titles1(i))
       daspect([1 1 1])
       %shading interp;
       view(2)
    end

    %nastavi pozicijo colormap legende
    hp4 = get(subplot(3,3,6),'Position');
    colorbar('Position', [1.33*hp4(1)  hp4(2)/2  hp4(3)/10  3*hp4(4)])

    
%     % narise odvisnost v eni smeri
    fig7=figure;
    set(fig7,'visible','off');
    for i = 1:9
       subplot(3,3,i);
       plot(eval(['lE', num2str(i)]));
       ylim([1.1*min(minar) 1.1*max(maxar)]);
       title(titles2(i))
    end

    %narise intenziteto
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
       %shading interp;
       view(2)
    end

    clear sE1 sE2 sE3 sE4 sE5 sE6 sE7 sE8 sE9
    clear lE1 lE2 lE3 lE4 lE5 lE6 lE7 lE8 lE9

    %WHERE TO SAVE THEM
    if isempty(file) == 0
        lam = abs(2*pi/value);
        out = strcat(file,'/',num2str(lam),'_',outgen,'_',whichplot,'_');
%         print(fig1,strcat(out,'H2D'),'-dpng');
%         pause(0.5)
%         print(fig2,strcat(out,'H2Dplane','.png'),'-dpng');
%         pause(1)
        print(fig3, strcat(out,'Hintensty', num2str(cycle), '.png'),'-dpng');
        pause(1) 
%         print(fig3, strcat(out,'H1D','.png'), '-dpng');
%         pause(0.5)        
%         print(fig4,strcat(out,'E2D_cycle',num2str(cycle),'.png'),'-dpng');
%         pause(1) 
%         print(fig5,strcat(out,'E2Dplane_n',num2str(cycle),'.png'),'-dpng');
%         pause(1) 
        print(fig6,strcat(out,'Eintensty',num2str(cycle),'.png'),'-dpng');
        pause(1) 
        print(fig7,strcat(out,'E1D_cycle',num2str(cycle),'.png'),'-dpng');
        pause(1) 
    end
end

clear fig1 fig2 fig3 fig4 fig5 fig6;

%SOME FUNCTIONS NEEDED IN CALCULATIONS
function f = ind(i,j,k)
    f = (k-1)*nx*ny + (j-1)*nx + (i-1) + 1;
end

function f = dHx(H,i,j,k)
    if i==1
        f = (H(2,j,k) - H(nx,j,k))/(2*dx);
    elseif i==nx
        f = (H(1,j,k) - H(nx-1,j,k))/(2*dx);
    else
        f = (H(i+1,j,k) - H(i-1,j,k))/(2*dx);
    end
end
function f = dHy(H,i,j,k)
     if j==1
        f = (H(i,2,k) - H(i,ny,k))/(2*dy);
    elseif j==ny
        f = (H(i,1,k) - H(i,ny-1,k))/(2*dy);
    else
        f = (H(i,j+1,k) - H(i,j-1,k))/(2*dy);
    end
end
function f = dHz(H,i,j,k)
    if k==1
        f = (H(i,j,2) - H(i,j,nz))/(2*dz);
    elseif k==nz
        f = (H(i,j,1) - H(i,j,nz-1))/(2*dz);
    else
        f = (H(i,j,k+1) - H(i,j,k-1))/(2*dz);
    end
end
function f = Imen(i,j,k)
    f = Exz(i,j,k)^2*Eyy(i,j,k) - 2*Exy(i,j,k)*Exz(i,j,k)*Eyz(i,j,k) ...
        + Exy(i,j,k)^2 * Ezz(i,j,k) + Exx(i,j,k)*(Eyz(i,j,k)^2 - ...
        Eyy(i,j,k)*Ezz(i,j,k));
end
%radius for droplets
function f = r(i,j,k)
       f = sqrt((i-nx/2)^2 + (j-ny/2)^2 + (k-nz/2)^2);
end


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

%radius for droplets
function f = r_qrt(i,j,k)
       f = sqrt((i)^2 + (j)^2 + (k)^2);
end

end

