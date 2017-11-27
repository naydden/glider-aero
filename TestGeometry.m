%Code to test the geometry
clc
clear all;

Nx=10;
Ny=50;
cr=1;
ct=0.5;
b=10;
m=0.02;
p=0.4;
sweep=10;
dihedral=5;
twist=15;
x_offset_W=0;
x_offset_T=5;
z_offset_W=0;
z_offset_T=0;

%Semi-Wing
[CoordLW,VortexLW,ControlPLW,DragPLW,NormalLW] = geometry (cr,ct,-b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset_W,z_offset_W);
[CoordRW,VortexRW,ControlPRW,DragPRW,NormalRW] = geometry (cr,ct,+b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset_W,z_offset_W);
%Semi-Tail
[CoordLT,VortexLT,ControlPLT,DragPLT,NormalLT] = geometry (0.5*cr,0.5*ct,-0.5*b,Nx,Ny,0,0,sweep,0,0,x_offset_T,z_offset_T);
[CoordRT,VortexRT,ControlPRT,DragPRT,NormalRT] = geometry (0.5*cr,0.5*ct,+0.5*b,Nx,Ny,0,0,sweep,0,0,x_offset_T,z_offset_T);
%Wing
CoordW = ensamblaje(CoordLW,CoordRW);
VortexW = ensamblaje(VortexLW,VortexRW);
ControlPW = ensamblaje(ControlPLW,ControlPRW);
DragPW = ensamblaje(DragPLW,DragPRW);
NormalW = ensamblaje(NormalLW,NormalRW);
%Tail
CoordT = ensamblaje(CoordLT,CoordRT);
VortexT = ensamblaje(VortexLT,VortexRT);
ControlPT = ensamblaje(ControlPLT,ControlPRT);
DragPT = ensamblaje(DragPLT,DragPRT);
NormalT = ensamblaje(NormalLT,NormalRT);
%Airplane
Coord = ensamblaje(CoordW,CoordT);
Vortex = ensamblaje(VortexW,VortexT);
ControlP = ensamblaje(ControlPW,ControlPT);
DragP = ensamblaje(DragPW,DragPT);
Normal = ensamblaje(NormalW,NormalT);

figure(1);
surf(CoordW(:,:,1),CoordW(:,:,2),CoordW(:,:,3));
hold on;
surf(CoordT(:,:,1),CoordT(:,:,2),CoordT(:,:,3));
axis equal;