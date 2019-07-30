function eps_chi_chi_mu = epsilon_bianisotropy_block(pdir,lambda2)
%THIS FUNCTION CALLS A FUNCTION TO BUILD EPSILON MATRIX AND BIANISOTTROPY MATRIX AND JOINS THEM IN
%A MATRIX i*[-EPSILON, BIANISOTROPY; BIANISOTROPY, MU]- 
global nx ny nz NO DN NOUT dirfield file R c DPMLs DPMLe gain;
[epsilon, mu] = epsilon_FB_block(pdir,lambda2);
bianisotropy = sparse(nx*ny*nz*3, nx*ny*nz*3);
eps_chi_chi_mu = i*[-epsilon, bianisotropy; bianisotropy, mu];

end
