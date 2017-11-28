function [Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr] = mirror (Coord,Vortex,ControlP,DragP,Normal)

Coord_Mirr=Coord;
Vortex_Mirr=Vortex;
ControlP_Mirr=ControlP;
DragP_Mirr=DragP;
Normal_Mirr=Normal;

Coord_Mirr(:,:,3)=-1*Coord_Mirr(:,:,3);
Vortex_Mirr(1:4,:,3)=-1*Vortex_Mirr(1:4,:,3);
ControlP_Mirr(:,3)=-1*ControlP_Mirr(:,3);
DragP_Mirr(:,3)=-1*DragP_Mirr(:,3);
Normal_Mirr(:,3)=-1*Normal_Mirr(:,3);

end