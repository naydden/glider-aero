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
twist=10;

[Coord,Vortex,ControlP,DragP,Normal] = total_geometry (cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist);

surf(Coord(:,:,1),Coord(:,:,2),Coord(:,:,3));
axis equal;