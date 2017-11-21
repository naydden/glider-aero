%Code to test the geometry
Nx=10;
Ny=50;
cr=1;
ct=0.5;
b=10;
m=0.02;
p=0.4;
flecha=20;
diedro=10;

x=zeros(Nx+1,Ny+1);
y=zeros(Nx+1,Ny+1);
z=zeros(Nx+1,Ny+1);
xp=zeros(Nx+1,Ny+1);
yp=zeros(Nx+1,Ny+1);
zp=zeros(Nx+1,Ny+1);
xc=zeros(Nx,Ny);
yc=zeros(Nx,Ny);
zc=zeros(Nx,Ny);

[x,y,z,xp,yp,zp,xc,yc,zc] = geometry (cr,ct,b,Nx,Ny,m,p,flecha,diedro);

surf(x,y,z);
axis equal;