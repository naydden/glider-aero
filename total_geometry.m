% It calculates the geometry of the whole wing
function [Coord,Vortex,ControlP,DragP,Normal] = total_geometry (cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist)

% Geometry of the left semi-wing
[Coord_left,CoordP_left,CoordC_left,CoordD_left,n_left] = geometry (cr,ct,-b/2,Nx,Ny/2,m,p,sweep,dihedral,twist);

% Geometry of the right semi-wing
[Coord_right,CoordP_right,CoordC_right,CoordD_right,n_right] = geometry (cr,ct,b/2,Nx,Ny/2,m,p,sweep,dihedral,twist);

N = Nx*Ny;
Coord=zeros(Nx+1,Ny+1,3);
Vortex=zeros(5,N,3);
ControlP=zeros(N,3);
DragP=zeros(N,3);
Normal=zeros(N,3);

% Assembly
for j = 1:Ny
    for i = 1:Nx
        if j<=Ny/2+1
            Coord(i,j,1) = Coord_left(i,j,1);
            Coord(i,j,2) = Coord_left(i,j,2);
            Coord(i,j,3) = Coord_left(i,j,3);
        else
            Coord(i,j,1) = Coord_right(i,j-Ny/2+1,1);
            Coord(i,j,2) = Coord_right(i,j-Ny/2+1,2);
            Coord(i,j,3) = Coord_right(i,j-Ny/2+1,3);
        end
    end
end
end