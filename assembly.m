function [Coord,Vortex,ControlP,DragP,Normal] = assembly(CoordA,VortexA,ControlPA,DragPA,NormalA,CoordB,VortexB,ControlPB,DragPB,NormalB)
Coord=[CoordA,CoordB];
Vortex=[VortexA,VortexB];
ControlP=[ControlPA;ControlPB];
DragP=[DragPA;DragPB];
Normal=[NormalA;NormalB];
end