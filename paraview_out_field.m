function [] = paraview_out_field(Vx,Vy,Vz,name)

%n ... which mode to draw
%pt --- transformation vector

global NX NY NZ out lambda

nx = NX;
ny = NY;
nz = NZ;

int = zeros(nx,ny,nz);

for i = 1:nx
    for j = 1:ny
        for k = 1:nz   
            int(i,j,k) = Vx(i,j,k)^2 + Vy(i,j,k)^2 + Vz(i,j,k)^2;
        end
    end
end

vecout = cat(4,Vx,Vy,Vz);
vecout = permute(vecout, [4,1,2,3]);

outvec = strcat(num2str(out),num2str(name),'_field_',num2str(nx),'_',num2str(ny),'_',num2str(nz),'_lda',num2str(lambda),'.raw');
outsq = strcat(num2str(out),num2str(name),'_int_',num2str(nx),'_',num2str(ny),'_',num2str(nz),'_lda',num2str(lambda),'.raw');

fid = fopen(outvec,'w');
fwrite(fid,vecout,'float');
fclose(fid);  

fid2 = fopen(outsq,'w');
fwrite(fid2,int,'float');
fclose(fid2);    

end


