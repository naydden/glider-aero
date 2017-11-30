clc;
clear all;
% Profile geometry
m_W = 0.02; p_W = 0.4;
% Wing geometry
cr_W = 1; ct_W = 1; b_W = 10; sweep_W = 0; dihedral_W = 0; twist_W = 0;
% Air
alpha = 3; x_offset_W = 0; z_offset_W=0; rho = 1.225;
Uinf = [1*cosd(alpha),0,1*sind(alpha)];
% Numerical
Nx = 2; Ny = 2;
deltaY = b_W/(2*Ny);
[Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);
% %Plane incidence
% incidence=3;
% [Coord,Vortex,ControlP,DragP,Normal] = rotation(Coord,Vortex,ControlP,DragP,Normal,incidence,0,cr_W,x_offset_W,z_offset_W);
% %Symetric plane
% [Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr] = mirror (Coord,Vortex,ControlP,DragP,Normal);
% %Ground-effect assembly
% [Coord,Vortex,ControlP,DragP,Normal] = assembly(Coord,Vortex,ControlP,DragP,Normal,Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr);

Gamma = circulation(Uinf,Vortex,ControlP,Normal);
[dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
[dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala');
L = lift(dLw,dLh,dLv);
M = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala');
Dind = drag(dDw,dDh,dDv);
[CL, CD, Cm] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M)