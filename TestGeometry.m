%Code to test the geometry
clc
clear all;

Nx=2;
Ny=4;
cr=2;
ct=2;
b=8;
m=0.02;
p=0.4;
sweep=0;
dihedral=0;
twist=0;

[Coord,Vortex,ControlP,DragP,Normal] = total_geometry (cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist);

surf(Coord(:,:,1),Coord(:,:,2),Coord(:,:,3));
axis equal;