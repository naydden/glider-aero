%Code to test the geometry
Nx=10;
Ny=50;
cr=1;
ct=0.5;
b=10;
m=0.02;
p=0.4;
sweep=20;
dihedral=10;

Coord=zeros(Nx+1,Ny+1,3);
CoordP=zeros(Nx+1,Ny+1,3);
CoordC=zeros(Nx,Ny,3);

[Coord,CoordP,CoordC] = geometry (cr,ct,b,Nx,Ny,m,p,sweep,dihedral);

surf(Coord(:,:,1),Coord(:,:,2),Coord(:,:,3));
axis equal;