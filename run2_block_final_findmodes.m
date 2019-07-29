lambda = 10;
lastmax = 0;

global file
file =  strcat('C:/Users/USER/Documents/Physics/fotonika');

%run multiple cycles of diagonalisation
for cycle=1:1

    global nx ny nz Lx Ly Lz dx dy dz NO DN NOUT nr dirfield pfield ref refz R DPMLe DPMLs oct c gain;

    %box size in px
    nx = 20; 
    ny = 20;
    nz = 20;
    %droplet radius
    R = 70;
    %box size in nm
    Lx = nx ;
    Ly = ny ;   %da vse merimo v enotah radija
    Lz = nz ;
    %pixel size
    dx = 1;
    dy = 1;
    dz = 1;
    % > 0 for gain inside droplet gain in pixels
    % already a combination of g * lambda
    gain = 0; %500 * 10^(-3);
    %ordinary refractive index and difference n_e - n_o
    NO = 1.2;%1.54;
    DN = 0.1;%0.17;
    NOUT = 1;%1.47;
    %thickness of PML, s = at start of domain (i=0) e = end of domain (i=n)
    %either or both can be applied
    DPMLs = [0,0,0];
    DPMLe = [0,0,0];
    %which octant is simulated - symmetry conditions for E
    oct = 5;
    %number of modes we want to find
    nr = 10;
    %YES if you want to plot director field
    pdir = 'YES';
    %YES for reflective BC
    ref = 'NO';
    refz = 'NO';
    BC1 = 'PEC';
    %YES if you want to plot electric and magnetic field
    pfield = 'YES';
    
    diagonalize= 'YES';
    %DRAWALL draws all the modes
    draw = 'DRAWALL';
    %dir field HELICONIC, HELICONICXY, HELICONICXZ, HELICONICYZ, RADIALD,
    %ISOTROPIC, ZERO ESCAPEDC FILE
    dirfield = 'ISOTROPIC';
    %centre of the droplet - ADD +1 FOR PML BOUNDARY
    c = [nx/2, ny/2, nz/2];

    disp(cycle)
    if lambda < 9
        break
    end

    sigma = (2*pi)^2/(lambda)^2;
    lambda =  2*pi/sqrt(sigma);
    
    outgen = strcat(num2str(nx),'x',num2str(ny),'x',num2str(nz));

    %create log file
    log = fopen(strcat(file,'/logfile_lda',num2str(lambda),'_',outgen,'.txt'),'w');
    logfile(log,lambda);  

    if strcmp(diagonalize, 'YES')    
        disp('Starting diagonalization process');  
        
        %write epsilon matrix
        B1 = epsilon_bianisotropy_block(pdir,lambda);
        
        %write derivative matrix
        A1 = FillA_S_block(BC1);
        
        %transformation to reduce bandwidth
        p = symrcm(A1);
        A = A1(p,p);
        B = B1(p,p);
        clear A1 B1 F M

        t1 = cputime;
        disp('Starting eigensolver');

        %EIGS solver - finds nr eigenvalues around sigma
        OPTS.tol = 1e-16;
        OPTS.maxit = 500;
        [V1,D] = eigs(A, B, nr, sigma);
        
        tdiag = cputime - t1;
        pt(p) = 1:length(p);
    end
        values = zeros(1,nr);
        vecH = zeros(3*nx*ny*nz,1);
        vecE = zeros(3*nx*ny*nz,1);
        
        fprintf(log,'\n Diagonalization complete\n');       

        %POSTPROCESS TOOLS
        for l = 1:nr
            values(1,l) = D(l,l);
            Q = 0;%real(sqrt(values(1,l))) / (2*imag(sqrt(values(1,l))));
            X = [num2str(l),'  ',num2str(values(1,l)),'  ',num2str(sqrt(4*pi^2/values(1,l))),'  ',num2str(Q)];
        end 

        for l = 1:nr     
            V = V1(:,l);
            V2 = V(pt);
            vecE(:,1) = V2(1:nx*ny*nz*3); 
            vecH(:,1) = V2(nx*ny*nz*3+1:nx*ny*nz*6);
            
            [~,~,~,~,dW] = script_func_block(vecH,vecE,values(1,l),'NOPML',l,outgen);

            lda = 2*pi/sqrt(values(1,l));
            Q = 0;%real(sqrt(values(1,l))) / (2*imag(sqrt(values(1,l))));
            
            %saves entire E eigenvector to a binary file
            if abs(Q) > 500
                outmode = strcat(file,'/Efield_lda',num2str(lda),'Q_',num2str(Q),'_',outgen,'.mat');
                save(outmode, 'vecE','-v7.3');
            end

            fprintf(log,'%d \t %f \t %f \t %f \t %f \n',l,values(1,l),lda,Q,dW);
        end 


        disp('Diagonalization complete')
        
    if max(values) > sigma
        sigma3 = max(values);
    else
        if max(values) > lastmax
            sigma3 = sigma;
        else
            sigma3 = sigma*1.05;
        end
    end
    sigma = sigma3;
    disp(sigma)
    if max(values)>lastmax
        lastmax = max(values);
    end
    
    lambda =  real(2*pi/sqrt(sigma)) * (1 - 1E-4);
    
    clearvars -except lambda lastmax sigma
        
end



function logfile(log,lambda)
    global nx ny nz Lx Ly Lz NO DN NOUT nr dirfield DPMLs DPMLe R c;
    fprintf(log,'Size (px): %d x %d x %d \n',nx,ny,nz);
    fprintf(log,'Box size (mum): %d x %d x %d \n',Lx,Ly,Lz);
    fprintf(log,'N_e: %3f \n',NO+DN);
    fprintf(log,'N_o: %3f \n',NO);
    fprintf(log,'N_out: %f \n \n',NOUT);  
    fprintf(log,'Nr. of modes: %d \n',nr);  
    fprintf(log,'Desired lambda (mum): %f \n \n',lambda);
    fprintf(log,'Director field: %s \n',dirfield);
    fprintf(log,'Radius: %f \n \n',R);
    fprintf(log,'DPML start: %f \t %f \t %f \n',DPMLs(1),DPMLs(2),DPMLs(3));
    fprintf(log,'DPML end: %f \t %f \t %f \n',DPMLe(1),DPMLe(2),DPMLe(3));    
    fprintf(log,'Drop center: %f \t %f \t %f \n',c(1),c(2),c(3));
    
end

%radius for droplets
function f = r_qrt(i,j,k)
       f = sqrt((i)^2 + (j)^2 + (k)^2);
end

