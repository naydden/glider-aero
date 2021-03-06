function [Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset,z_offset)

[CoordLW,VortexLW,ControlPLW,DragPLW,NormalLW] = geometry (cr,ct,-b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset,z_offset);
[CoordRW,VortexRW,ControlPRW,DragPRW,NormalRW] = geometry (cr,ct,+b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset,z_offset);

[Coord,Vortex,ControlP,DragP,Normal] = assembly(CoordLW,VortexLW,ControlPLW,DragPLW,NormalLW,CoordRW,VortexRW,ControlPRW,DragPRW,NormalRW);

end