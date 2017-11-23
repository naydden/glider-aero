%Code to test the geometry
clc
clear all;

Nx=10;
Ny=50;
cr=1;
ct=1;
b=-10;
m=0.02;
p=0.4;
sweep=20;
dihedral=10;
twist=10;

Coord=zeros(Nx+1,Ny+1,3);
CoordP=zeros(Nx+1,Ny+1,3);
CoordC=zeros(Nx,Ny,3);
n=zeros(Nx,Ny,3);

[Coord,CoordP,CoordC,CoordD,n] = geometry (cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist);

surf(Coord(:,:,1),Coord(:,:,2),Coord(:,:,3));
axis equal;