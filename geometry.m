function [Coord,CoordP,CoordC,CoordD,n] = geometry (cr,ct,b,Nx,Ny,m,p,sweep,dihedral)

x=zeros(Nx+1,Ny+1);
y=zeros(Nx+1,Ny+1);
z=zeros(Nx+1,Ny+1);
xp=zeros(Nx+1,Ny+1);
yp=zeros(Nx+1,Ny+1);
zp=zeros(Nx+1,Ny+1);
xc=zeros(Nx,Ny);
yc=zeros(Nx,Ny);
zc=zeros(Nx,Ny);
n=zeros(Nx,Ny,3);

for i=1:Nx+1
    y(i,:)=linspace(0,b/2,Ny+1); %Y Coordinate of the points of the elemnts
    yp(i,:)=linspace(0,b/2,Ny+1); %Y Coordinate of the points of the vortex
end

for j=1:Ny+1
    c=cr-(cr-ct)*y(1,j)/(b/2); %Computation of the chord of the wing segment
    for i=1:Nx+1
        x(i,j)=cr/4-c/4+c*(i-1)/Nx+y(i,j)*tand(sweep); %X Coordinate of the points of the elemnts
        xp(i,j)=cr/4-c/4+c*(i-1)/Nx+c/(4*Nx)+yp(i,j)*tand(sweep); %X Coordinate of the points of the vortex
        xadim=(x(i,j)-x(1,j))/c; %Adimensionalization of the x coordinatie of the points
        z(i,j)=camber(xadim,m,p)+y(i,j)*tand(dihedral); %Z Coordinate of the points of the elemnts
        xadimp=(xp(i,j)-x(1,j))/c; %Adimensionalization of the x coordinatie of the vector
        zp(i,j)=camber(xadim,m,p)+yp(i,j)*tand(dihedral); %Z Coordinate of the points of the vortex
        if c==0
            z(i,j)=y(i,j)*tand(dihedral);
            zp(i,j)=yp(i,j)*tand(dihedral);
        end
    end
end

for j=1:Ny
    for i=1:Nx
        yc(i,j)=(y(i,j)+y(i,j+1))/2; %Y Coordinate  of the control points
    end
    c=cr-(cr-ct)*yc(1,j)/(b/2); %Chord of the wing segment
    for i=1:Nx
        xc(i,j)=cr/4-c/4+c*(i-1)/Nx+c*3/(4*Nx)+yc(i,j)*tand(sweep); %X Coordinate  of the control points
        xd(i,j)=(xp(i,j)+xp(i,j+1))/2;  % X Coordinate of the drag control points
        xadimc=(xc(i,j)-((x(1,j)+x(1,j+1))/2))/c; %Adimensionalization  of the x coordinate of the control point
        xadimd=(xd(i,j)-((x(1,j)+x(1,j+1))/2))/c; %Adimensionalization  of the x coordinate of the drag control point
        zc(i,j)=camber(xadimc,m,p)+yc(i,j)*tand(dihedral); %Z Coordinate  of the control points
        zd(i,j)=camber(xadimd,m,p)+yc(i,j)*tand(dihedral); %Z Coordinate of the drag control points
    end
end

for i=1:Nx
    for j=1:Ny
        AC=[xp(i+1,j+1)-xp(i,j) yp(i+1,j+1)-yp(i,j) zp(i+1,j+1)-zp(i,j)];
        DB=[xp(i,j+1)-xp(i+1,j) yp(i,j+1)-yp(i+1,j) zp(i,j+1)-zp(i+1,j)];
        n(i,j,:)=cross(AC,DB)/norm(cross(AC,DB));
    end
end

    Coord(:,:,1)=x(:,:);
    Coord(:,:,2)=y(:,:);
    Coord(:,:,3)=z(:,:);
    CoordP(:,:,1)=xp(:,:);
    CoordP(:,:,2)=yp(:,:);
    CoordP(:,:,3)=zp(:,:);
    CoordC(:,:,1)=xc(:,:);
    CoordC(:,:,2)=yc(:,:);
    CoordC(:,:,3)=zc(:,:);
    CoordD(:,:,1)=xd(:,:);
    CoordD(:,:,2)=yc(:,:);
    CoordD(:,:,3)=zd(:,:);
    
end