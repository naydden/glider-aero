function [c,x,y] = geometry (cr,ct,b,Nx,Ny)
for j=1:Ny
    c(j)=cr-(cr-ct)/Ny*(j-1);
    for i=1:Nx
        x(i,j)=cr/4-c(j)/4+c(j)/Nx*(i-1);
        y(i,j)=b/(2*Ny)*(j-1);
    end
end
end