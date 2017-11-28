%Code to test the geometry
clc
clear all;

Nx=2;
Ny=4;
cr=3;
ct=1;
b=10;
m=0.02;
p=0.4;
sweep=20;
dihedral=40;
twist=50;
x_offset=0;
z_offset=2;
incidence=20;
lateral=0;

[CoordW,VortexW,ControlPW,DragPW,NormalW] = wing_assembly (cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset,z_offset);
[CoordH,VortexH,ControlPH,DragPH,NormalH] = wing_assembly (2,0,0.5*b,Nx,Ny,0,0,50,10,0,3,z_offset);
[CoordV,VortexV,ControlPV,DragPV,NormalV] = geometry (2,1,0.1*b,Nx,Ny,0,0,60,0,0,3,z_offset);
[CoordV,VortexV,ControlPV,DragPV,NormalV] = rotation(CoordV,VortexV,ControlPV,DragPV,NormalV,0,90,2,3,z_offset);

[CoordT,VortexT,ControlPT,DragPT,NormalT] = assembly(CoordH,VortexH,ControlPH,DragPH,NormalH,CoordV,VortexV,ControlPV,DragPV,NormalV);
[Coord,Vortex,ControlP,DragP,Normal] = assembly(CoordW,VortexW,ControlPW,DragPW,NormalW,CoordT,VortexT,ControlPT,DragPT,NormalT);

[Coord,Vortex,ControlP,DragP,Normal] = rotation(Coord,Vortex,ControlP,DragP,Normal,incidence,lateral,cr,x_offset,z_offset);

[Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr] = mirror (Coord,Vortex,ControlP,DragP,Normal);

figure(1);
surf(Coord(:,1:2*(Ny+1),1),Coord(:,1:2*(Ny+1),2),Coord(:,1:2*(Ny+1),3));
hold on;
surf(Coord(:,2*(Ny+1)+1:4*(Ny+1),1),Coord(:,2*(Ny+1)+1:4*(Ny+1),2),Coord(:,2*(Ny+1)+1:4*(Ny+1),3));
hold on;
surf(Coord(:,4*(Ny+1)+1:5*(Ny+1),1),Coord(:,4*(Ny+1)+1:5*(Ny+1),2),Coord(:,4*(Ny+1)+1:5*(Ny+1),3));
hold on;
surf(Coord_Mirr(:,1:2*(Ny+1),1),Coord_Mirr(:,1:2*(Ny+1),2),Coord_Mirr(:,1:2*(Ny+1),3));
hold on;
surf(Coord_Mirr(:,2*(Ny+1)+1:4*(Ny+1),1),Coord_Mirr(:,2*(Ny+1)+1:4*(Ny+1),2),Coord_Mirr(:,2*(Ny+1)+1:4*(Ny+1),3));
hold on;
surf(Coord_Mirr(:,4*(Ny+1)+1:5*(Ny+1),1),Coord_Mirr(:,4*(Ny+1)+1:5*(Ny+1),2),Coord_Mirr(:,4*(Ny+1)+1:5*(Ny+1),3));
axis equal;