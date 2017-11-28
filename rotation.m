function [CoordR,VortexR,ControlPR,DragPR,NormalR] = rotation(Coord,Vortex,ControlP,DragP,Normal,incidence,lateral,cr,x_offset,z_offset)

CoordR=Coord;
VortexR=Vortex;
ControlPR=ControlP;
DragPR=DragP;
NormalR=Normal;

%Set the origin at x=cr/4, z=0
CoordR(:,:,1)=Coord(:,:,1)-(x_offset+0.25*cr);
VortexR(1:4,:,1)=Vortex(1:4,:,1)-(x_offset+0.25*cr);
ControlPR(:,1)=ControlP(:,1)-(x_offset+0.25*cr);
DragPR(:,1)=DragP(:,1)-(x_offset+0.25*cr);

CoordR(:,:,3)=Coord(:,:,3)-z_offset;
VortexR(1:4,:,3)=Vortex(1:4,:,3)-z_offset;
ControlPR(:,3)=ControlP(:,3)-z_offset;
DragPR(:,3)=DragP(:,3)-z_offset;

%Incidence rotation
for i=1:size(Coord,1)
    for j=1:size(Coord,2)
        Coord(i,j,1)=CoordR(i,j,1)*cosd(incidence)+CoordR(i,j,3)*sind(incidence);
        Coord(i,j,3)=CoordR(i,j,3)*cosd(incidence)-CoordR(i,j,1)*sind(incidence);
    end
end

for i=1:size(ControlP,1)
    Vortex(1:4,i,1)=VortexR(1:4,i,1)*cosd(incidence)+VortexR(1:4,i,3)*sind(incidence);
    Vortex(1:4,i,3)=VortexR(1:4,i,3)*cosd(incidence)-VortexR(1:4,i,1)*sind(incidence);
    ControlP(i,1)=ControlPR(i,1)*cosd(incidence)+ControlPR(i,3)*sind(incidence);
    ControlP(i,3)=ControlPR(i,3)*cosd(incidence)-ControlPR(i,1)*sind(incidence);
    DragP(i,1)=DragPR(i,1)*cosd(incidence)+DragPR(i,3)*sind(incidence);
    DragP(i,3)=DragPR(i,3)*cosd(incidence)-DragPR(i,1)*sind(incidence);
    Normal(i,1)=NormalR(i,1)*cosd(incidence)+NormalR(i,3)*sind(incidence);
    Normal(i,3)=NormalR(i,3)*cosd(incidence)-NormalR(i,1)*sind(incidence);
end

%Lateral rotation
for i=1:size(Coord,1)
    for j=1:size(Coord,2)
        CoordR(i,j,2)=Coord(i,j,2)*cosd(lateral)-Coord(i,j,3)*sind(lateral);
        CoordR(i,j,3)=Coord(i,j,3)*cosd(lateral)+Coord(i,j,2)*sind(lateral);
    end
end

for i=1:size(ControlP,1)
    VortexR(1:4,i,2)=Vortex(1:4,i,2)*cosd(lateral)+Vortex(1:4,i,3)*sind(lateral);
    VortexR(1:4,i,3)=Vortex(1:4,i,3)*cosd(lateral)-Vortex(1:4,i,2)*sind(lateral);
    ControlPR(i,2)=ControlP(i,2)*cosd(lateral)+ControlP(i,3)*sind(lateral);
    ControlPR(i,3)=ControlP(i,3)*cosd(lateral)-ControlP(i,2)*sind(lateral);
    DragPR(i,2)=DragP(i,2)*cosd(lateral)+DragP(i,3)*sind(lateral);
    DragPR(i,3)=DragP(i,3)*cosd(lateral)-DragP(i,2)*sind(lateral);
    NormalR(i,2)=Normal(i,2)*cosd(lateral)+Normal(i,3)*sind(lateral);
    NormalR(i,3)=Normal(i,3)*cosd(lateral)-Normal(i,2)*sind(lateral);
end

%Repositio the wing at original coordinates
CoordR(:,:,1)=CoordR(:,:,1)+(x_offset+0.25*cr);
VortexR(1:4,:,1)=VortexR(1:4,:,1)+(x_offset+0.25*cr);
ControlPR(:,1)=ControlPR(:,1)+(x_offset+0.25*cr);
DragPR(:,1)=DragPR(:,1)+(x_offset+0.25*cr);

CoordR(:,:,3)=CoordR(:,:,3)+z_offset;
VortexR(1:4,:,3)=VortexR(1:4,:,3)+z_offset;
ControlPR(:,3)=ControlPR(:,3)+z_offset;
DragPR(:,3)=DragPR(:,3)+z_offset;
end