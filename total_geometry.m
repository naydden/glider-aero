% It calculates the geometry of the whole wing
function [Coord,Vortex,ControlC,DragP,Normal] = total_geometry (cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist)
% b: wingspan of the WHOLE wing
% Ny: number of elements of the WHOLE wing (it has to be an even number)

% Geometry of the left semi-wing
[Coord_left,Vortex_left,ControlC_left,DragP_left,Normal_left] = geometry (cr,ct,-b,Nx,Ny/2,m,p,sweep,dihedral,twist);

% Geometry of the right semi-wing
[Coord_right,Vortex_right,ControlC_right,DragP_right,Normal_right] = geometry (cr,ct,b,Nx,Ny/2,m,p,sweep,dihedral,twist);

N = Nx*Ny;
Coord=zeros(Nx+1,Ny+1,3);
Vortex=zeros(5,N,3);
ControlC=zeros(N,3);
DragP=zeros(N,3);
Normal=zeros(N,3);

% Assembly
for j = 1:Ny+1
    for i = 1:Nx+1
        if j<=Ny/2+1
            Coord(i,j,:) = Coord_left(i,j,:);
        else
            Coord(i,j,:) = Coord_right(i,j-Ny/2,:);
        end
    end
end

for i = 1:N
    if i<=N/2
        ControlC(i,:) = ControlC_left(i,:);
        DragP(i,:) = DragP_left(i,:);
        Normal(i,:) = Normal_left(i,:);
        Vortex(:,i,:) = Vortex_left(:,i,:);
    else
        ControlC(i,:) = ControlC_right(i-N/2,:);
        DragP(i,:) = DragP_right(i-N/2,:);
        Normal(i,:) = Normal_right(i-N/2,:);
        Vortex(:,i,:) = Vortex_right(:,i-N/2,:);
    end
end