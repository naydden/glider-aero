function [x,y,xp,yp,xc,yc] = geometry (cr,ct,b,Nx,Ny)

for i=1:Nx+1
    y(i,:)=linspace(0,b/2,Ny+1);
    yp(i,:)=linspace(0,b/2,Ny+1);
end

for j=1:Ny+1
    c=cr-(cr-ct)*y(1,j)/(b/2);
    for i=1:Nx+1
        x(i,j)=cr/4-c/4+c*(i-1)/Nx;
        xp(i,j)=cr/4-c/4+c*(i-1)/Nx+c/(4*Nx);
    end
end

for j=1:Ny
    yc(:,j)=(y(i,j)+y(i,j+1))/2;
    c=cr-(cr-ct)*yc(1,j)/(b/2);
    for i=1:Nx
        xc(i,j)=cr/4-c/4+c*(i-1)/Nx+c*3/(4*Nx);
    end
end

end