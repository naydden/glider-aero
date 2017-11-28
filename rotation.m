function [CoordR,VortexR,ControlPR,DragPR,NormalR] = rotation(Coord,Vortex,ControlP,DragP,Normal,alpha,cr,x_offset,z_offset)

CoordR=Coord;
VortexR=Vortex;
ControlPR=ControlP;
DragPR=DragP;
NormalR=Normal;

Coord(:,:,1)=Coord(:,:,1)-(x_offset+0.25*cr);
Vortex(1:4,:,1)=Vortex(1:4,:,1)-(x_offset+0.25*cr);
ControlP(:,1)=ControlP(:,1)-(x_offset+0.25*cr);
DragP(:,1)=DragP(:,1)-(x_offset+0.25*cr);

Coord(:,:,3)=Coord(:,:,3)-z_offset;
Vortex(1:4,:,3)=Vortex(1:4,:,3)-z_offset;
ControlP(:,3)=ControlP(:,3)-z_offset;
DragP(:,3)=DragP(:,3)-z_offset;

for i=1:size(Coord,1)
    for j=1:size(Coord,2)
        CoordR(i,j,1)=Coord(i,j,1)*cosd(alpha)+Coord(i,j,3)*sind(alpha);
        CoordR(i,j,3)=Coord(i,j,3)*cosd(alpha)-Coord(i,j,1)*sind(alpha);
    end
end

for i=1:size(ControlP,1)
    Vortex(1:4,i,1)=Vortex(1:4,i,1)*cosd(alpha)+Vortex(1:4,i,3)*sind(alpha);
    Vortex(1:4,i,3)=Vortex(1:4,i,3)*cosd(alpha)-Vortex(1:4,i,1)*sind(alpha);
    ControlP(i,1)=ControlP(i,1)*cosd(alpha)+ControlP(i,3)*sind(alpha);
    ControlP(i,3)=ControlP(i,3)*cosd(alpha)-ControlP(i,1)*sind(alpha);
    DragP(i,1)=DragP(i,1)*cosd(alpha)+DragP(i,3)*sind(alpha);
    DragP(i,3)=DragP(i,3)*cosd(alpha)-DragP(i,1)*sind(alpha);
    Normal(i,1)=Normal(i,1)*cosd(alpha)+Normal(i,3)*sind(alpha);
    Normal(i,3)=Normal(i,3)*cosd(alpha)-Normal(i,1)*sind(alpha);
end

CoordR(:,:,1)=CoordR(:,:,1)+(x_offset+0.25*cr);
VortexR(1:4,:,1)=Vortex(1:4,:,1)+(x_offset+0.25*cr);
ControlPR(:,1)=ControlP(:,1)+(x_offset+0.25*cr);
DragPR(:,1)=DragP(:,1)+(x_offset+0.25*cr);

CoordR(:,:,3)=CoordR(:,:,3)+z_offset;
VortexR(1:4,:,3)=Vortex(1:4,:,3)+z_offset;
ControlPR(:,3)=ControlP(:,3)+z_offset;
DragPR(:,3)=DragP(:,3)+z_offset;
end