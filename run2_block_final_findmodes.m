lambda = 100;
lastmax = 0;

global file
%file =  strcat('/data/Matlab/bianisotropyFB/Alrods_air/kvecG');
file =  strcat('C:/Users/USER/Documents/Physics/fotonika');
mkdir(strcat(file,'/dir'));

kx = [0 1/4 2/4 3/4 1 1   1   1   1 3/4 2/4 1/4 0];
ky = [0 0   0   0   0 1/4 2/4 3/4 1 3/4 2/4 1/4 0];

omegakTE = fopen(strcat(file,'/omegakTE.txt'),'w');
omegakTM = fopen(strcat(file,'/omegakTM.txt'),'w');



%run multiple cycles of diagonalisation
for cycle=1:1

    global nx ny nz u1 u2 u3 Lx Ly Lz dx dy dz NO DN NOUT nr pfield ref refz R DPMLe DPMLs oct c gain Delta modetype kvec;

    nx = 15; %box size in pixels
    ny = 15;
    nz = 10;
    
    %triangular grid
    u1 = [1, 0, 0];
    u2 = [0.5, 0.5*sqrt(3), 0];
    u3 = [0, 0, 1];
    
    R = nx/5; %droplet radius
    Lx = nx ;
    Ly = ny ;
    Lz = nz ;
    %pixel size
    dx = 1;
    dy = 1;
    dz = 1;
    % > 0 for gain inside droplet gain in pixels
    % already a combination of g * lambda
    gain = 0; %500 * 10^(-3);
    Delta = 0.1; %parameter in Khanikaev bianisotropy
    %ordinary refractive index and difference n_e - n_o
    NO = sqrt(8.9); %1.54;
    DN = 0.0; %0.17;
    NOUT = 1; %1.47;
    %thickness of PML, s = at start of domain (i=0) e = end of domain (i=n)
    %either or both can be applied
    DPMLs = [0,0,0];
    DPMLe = [0,0,0];
    %which octant is simulated - symmetry conditions for E
    oct = 5;
    %number of modes we want to find
    nr = 10;
    %YES for reflective BC
    ref = 'NO';
    refz = 'NO';
    BC1 = 'PMC';
    %YES if you want to plot electric and magnetic field
    pfield = 'YES';

    diagonalize = 'YES';
    %DRAWALL draws all the modes
    draw = 'DRAWALL';
    
    epstype = 'CYLINDER'; %dir field HELICONIC, HELICONICXY, HELICONICXZ, HELICONICYZ, RADIALD,
    %                                  ISOTROPIC, ZERO, ESCAPEDC, FILE, KHANIKAEV
    mutype = 'IDENTITY';
    biantype = 'ZERO';
    
    modetype = 'ALL';
    
    %centre of the droplet - ADD +1 FOR PML BOUNDARY
    c = [nx*0.75, (ny/nx)*nx*sqrt(3)/4, nz/2];
    %wave vector in units omega/2*pi*c
    kvec = pi*[kx(cycle) ky(cycle) 0]; 

    sigma = 2*pi/lambda;

    outgen = strcat(num2str(nx),'x',num2str(ny),'x',num2str(nz));
    
    log = fopen(strcat(file,'/logfile_lda',num2str(lambda),'_',outgen,'.txt'),'w');
    logfile(log,lambda);

    if strcmp(diagonalize, 'YES')    
        disp('Starting diagonalization process');  

        %write epsilon matrix
        B1 = epsilon_bianisotropy_mu_block(epstype, biantype, mutype, lambda);

        %write derivative matrix
        A1 = FillA_FB_block();

        %transformation to reduce bandwidth
        %p = symrcm(A1);
        A = A1;%(p,p);
        B = B1;%(p,p);
        clear A1 B1 F M

        t1 = cputime;
        disp('Starting eigensolver');

        %EIGS solver - finds nr eigenvalues around sigma
        OPTS.tol = 1e-16;
        OPTS.maxit = 1000;
        [V1,D] = eigs(A,B, nr, sigma, OPTS);

        tdiag = cputime - t1;
        %pt(p) = 1:length(p);
    end
        values = zeros(1,nr);
        vecH = zeros(3*nx*ny*nz,1);
        vecE = zeros(3*nx*ny*nz,1);

        %POSTPROCESS TOOLS
        for l = 1:nr
             values(1,l) = D(l,l);
        end 
        
        for l = 1:nr     
            V = V1(:,l);
            V2 = V;%(pt);
            vecE(:,1) = V2(1:nx*ny*nz*3); 
            vecH(:,1) = V2(nx*ny*nz*3+1:nx*ny*nz*6);

            [~,~,~,~,dW,mH,mE,type] = script_func_block(vecH,vecE,values(1,l),'NOPML',l,outgen);

            lda = 2*pi/values(1,l)
            
            fprintf(log,'%d \t %f \t %f \t %s \t %f \t %f \t %f \t %f \t %f \t %f \n',l,values(1,l)*nx/(2*pi),lda,type,mH(1),mH(2),mH(3),mE(1),mE(2),mE(3));
            if strcmp(type, 'TE')   
                fprintf(omegakTE,'%d \t %f \t %f \t %s \t %f \t %f \t %f \t %f \t %f \t %f \n',cycle,values(1,l)*nx/(2*pi),lda,type,mH(1),mH(2),mH(3),mE(1),mE(2),mE(3));
            elseif strcmp(type, 'TM')   
                fprintf(omegakTM,'%d \t %f \t %f \t %s \t %f \t %f \t %f \t %f \t %f \t %f \n',cycle,values(1,l)*nx/(2*pi),lda,type,mH(1),mH(2),mH(3),mE(1),mE(2),mE(3));
            else
                continue
            end
        end 


        disp('Diagonalization complete')

%     if max(values) > sigma
%         sigma3 = max(values);
%     else
%         if max(values) > lastmax
%             sigma3 = sigma;
%         else
%             sigma3 = sigma*1.05;
%         end
%     end
%     sigma = sigma3;
%     disp(sigma)
%     if max(values)>lastmax
%         lastmax = max(values);
%     end
% 
%     lambda =  real(2*pi/sigma) * (1 - 1E-4);
% 
     clearvars -except lambda sigma ky kx omegakTE omegakTM file

end

function logfile(log,lambda)
    global nx ny nz Lx Ly Lz NO DN NOUT nr dirfield DPMLs DPMLe R c kvec;
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
    fprintf(log,'Drop center: %f \t %f \t %f \n \n',c(1),c(2),c(3));
    fprintf(log,'k vector: %f \t %f \t %f \n',kvec(1),kvec(2),kvec(3));
    
end

