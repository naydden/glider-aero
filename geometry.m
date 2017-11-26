%Function that computs the geometric points and vectors needed for the
%program
function [Coord,Vortex,ControlP,DragP,Normal] = geometry (cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist)

x=zeros(Nx+1,Ny+1);
y=zeros(Nx+1,Ny+1);
z=zeros(Nx+1,Ny+1);
xp=zeros(Nx+1,Ny+1);
yp=zeros(Nx+1,Ny+1);
zp=zeros(Nx+1,Ny+1);
xc=zeros(Nx,Ny);
yc=zeros(Nx,Ny);
zc=zeros(Nx,Ny);
xd=zeros(Nx,Ny);
zd=zeros(Nx,Ny);
n=zeros(Nx,Ny,3);

Coord=zeros(Nx+1,Ny+1,3);
Vortex=zeros(5,Nx*Ny,3);
ControlP=zeros(Nx*Ny,3);
DragP=zeros(Nx*Ny,3);
Normal=zeros(Nx*Ny,3);

%Geometry points and vortex points
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
        zp(i,j)=camber(xadimp,m,p)+yp(i,j)*tand(dihedral); %Z Coordinate of the points of the vortex
        if c==0
            z(i,j)=y(i,j)*tand(dihedral);
            zp(i,j)=yp(i,j)*tand(dihedral);
        end
    end
end

%Control points for lift and drag
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

%Twist for geomety and vortex points
for i=1:Nx+1
    for j=1:Ny+1
        twisti=twist*y(i,j)/(b/2);
        xtwist=0.25*cr+(x(i,j)-0.25*cr)*cosd(twisti)+z(i,j)*sind(twisti);
        ztwist=z(i,j)*cosd(twisti)-(x(i,j)-0.25*cr)*sind(twisti);
        x(i,j)=xtwist;
        z(i,j)=ztwist;
        xtwist=0.25*cr+(xp(i,j)-0.25*cr)*cosd(twisti)+zp(i,j)*sind(twisti);
        ztwist=zp(i,j)*cosd(twisti)-(xp(i,j)-0.25*cr)*sind(twisti);
        xp(i,j)=xtwist;
        zp(i,j)=ztwist;
    end
end

%Twist for control points
for i=1:Nx
    for j=1:Ny
        twisti=twist*yc(i,j)/(b/2);
        xtwist=0.25*cr+(xc(i,j)-0.25*cr)*cosd(twisti)+zc(i,j)*sind(twisti);
        ztwist=zc(i,j)*cosd(twisti)-(xc(i,j)-0.25*cr)*sind(twisti);
        xc(i,j)=xtwist;
        zc(i,j)=ztwist;
        xtwist=0.25*cr+(xd(i,j)-0.25*cr)*cosd(twisti)+zd(i,j)*sind(twisti);
        ztwist=zd(i,j)*cosd(twisti)-(xd(i,j)-0.25*cr)*sind(twisti);
        xd(i,j)=xtwist;
        zd(i,j)=ztwist;
    end
end

%Normal vector
for i=1:Nx
    for j=1:Ny
        AC=[xp(i+1,j+1)-xp(i,j) yp(i+1,j+1)-yp(i,j) zp(i+1,j+1)-zp(i,j)];
        DB=[xp(i,j+1)-xp(i+1,j) yp(i,j+1)-yp(i+1,j) zp(i,j+1)-zp(i+1,j)];
        ni=cross(AC,DB);
        ni=ni/norm(ni);
        n(i,j,:)=ni;
    end
end

Coord(:,:,1)=x(:,:);
Coord(:,:,2)=y(:,:);
Coord(:,:,3)=z(:,:);

for i=1:Nx
    for j=1:Ny
        Vortex(1,(i-1)*Ny+j,:)=[xp(i,j) yp(i,j) zp(i,j)];
        Vortex(2,(i-1)*Ny+j,:)=[xp(i,j+1) yp(i,j+1) zp(i,j+1)];
        Vortex(3,(i-1)*Ny+j,:)=[xp(i+1,j+1) yp(i+1,j+1) zp(i+1,j+1)];
        Vortex(4,(i-1)*Ny+j,:)=[xp(i+1,j) yp(i+1,j) zp(i+1,j)];
        ControlP((i-1)*Ny+j,:)=[xc(i,j) yc(i,j) zc(i,j)];
        DragP((i-1)*Ny+j,:)=[xd(i,j) yc(i,j) zd(i,j)];
        Normal((i-1)*Ny+j,:)=n(i,j,:);
        if i==Nx
            Vortex(5,(i-1)*Ny+j,:)=1;
        else
            Vortex(5,(i-1)*Ny+j,:)=0;
        end
    end
end
    
end