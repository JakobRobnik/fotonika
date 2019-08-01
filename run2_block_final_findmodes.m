lambda = 30;
lastmax = 0;

global file
file =  strcat('C:/Users/USER/Documents/Physics/fotonika');

%run multiple cycles of diagonalisation
for cycle=1:1

    global nx ny nz Lx Ly Lz dx dy dz NO DN NOUT nr pfield ref refz R DPMLe DPMLs oct c gain;

    nx = 21; %box size in pixels
    ny = 21;
    nz = 5;
    
    R = 70; %droplet radius
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
    %ordinary refractive index and difference n_e - n_o
    NO = 1.2; %1.54;
    DN = 0.3; %0.17;
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

    diagonalize = 'YES'
    %DRAWALL draws all the modes
    draw = 'DRAWALL';
    
    epstype = 'HELICONIC'; %dir field HELICONIC, HELICONICXY, HELICONICXZ, HELICONICYZ, RADIALD,
    %                                  ISOTROPIC, ZERO, ESCAPEDC, FILE, KHANIKAEV
    mutype = 'HELICONIC';
    biantype = 'KHANIKAEV';
    
    %centre of the droplet - ADD +1 FOR PML BOUNDARY
    c = [nx/2, ny/2, nz/2];


    sigma = 2*pi/lambda;

    outgen = strcat(num2str(nx),'x',num2str(ny),'x',num2str(nz));

    if strcmp(diagonalize, 'YES')    
        disp('Starting diagonalization process');  

        %write epsilon matrix
        B1 = epsilon_bianisotropy_block(pdir, epstype, biantype, mutype, lambda);

        %write derivative matrix
        A1 = FillA_FB_block(BC1);

        %transformation to reduce bandwidth
        %p = symrcm(A1);
        A = A1;%(p,p);
        B = B1;%(p,p);
        clear A1 B1 F M

        t1 = cputime;
        disp('Starting eigensolver');

        %EIGS solver - finds nr eigenvalues around sigma
        OPTS.tol = 1e-16;
        OPTS.maxit = 500;
        [V1,D] = eigs(A, B, nr, sigma);

        tdiag = cputime - t1;
        %pt(p) = 1:length(p);
    end
        values = zeros(1,nr);
        vecH = zeros(3*nx*ny*nz,1);
        vecE = zeros(3*nx*ny*nz,1);

        %POSTPROCESS TOOLS
        for l = 1:nr
             values(1,l) = D(l,l);
            Q = 0;%real(sqrt(values(1,l))) / (2*imag(sqrt(values(1,l))));
            X = [num2str(l),'  ',num2str(values(1,l)),'  ',num2str(2*pi/values(1,l)),'  ',num2str(Q)];
        end 
        
        for l = 1:nr     
            V = V1(:,l);
            V2 = V;%(pt);
            vecE(:,1) = V2(1:nx*ny*nz*3); 
            vecH(:,1) = V2(nx*ny*nz*3+1:nx*ny*nz*6);

            [~,~,~,~,dW] = script_func_block(vecH,vecE,values(1,l),'NOPML',l,outgen);

            lda = 2*pi/values(1,l)
            Q = 0;%real(sqrt(values(1,l))) / (2*imag(sqrt(values(1,l))));

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

    lambda =  real(2*pi/sigma) * (1 - 1E-4);

    clearvars -except lambda lastmax sigma

end


%radius for droplets
function f = r_qrt(i,j,k)
       f = sqrt((i)^2 + (j)^2 + (k)^2);
end

%martina lebar
