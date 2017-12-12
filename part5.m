% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0;
Uinf = [1,0,0];

%Wing
cr_W=1; ct_W=1*cr_W; b_W=10;
sweep_W=0; dihedral_W=0; twist_W=0;
x_offset_W=0; z_offset_W=1;

[CoordW,VortexW,ControlPW,DragPW,NormalW] = wing_assembly (cr_W,ct_W,b_W,Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);

%Horizontal tail
cr_H=0.5*cr_W; ct_H=1*cr_H; b_H=0.25*b_W;
sweep_H=0; dihedral_H=0; twist_H=0;
x_offset_H=4+0.25*cr_W-0.25*cr_H; z_offset_H=z_offset_W;

[CoordH,VortexH,ControlPH,DragPH,NormalH] = wing_assembly (cr_H,ct_H,b_H,Nx,Ny,m_H,p_H,sweep_H,dihedral_H,twist_H,x_offset_H,z_offset_H);

%Vertical tail
cr_V=1*cr_H; ct_V=1*cr_V; b_V=2/3*b_H;
sweep_V=0; dihedral_V=0; twist_V=0;
x_offset_V=x_offset_H; z_offset_V=z_offset_W;

[CoordV,VortexV,ControlPV,DragPV,NormalV] = geometry (cr_V,ct_V,b_V,Nx,Ny,m_V,p_V,sweep_V,dihedral_V,twist_V,x_offset_V,z_offset_V);
[CoordV,VortexV,ControlPV,DragPV,NormalV] = rotation(CoordV,VortexV,ControlPV,DragPV,NormalV,0,90,cr_V,x_offset_V,z_offset_V);

%Tail assembly
[CoordT,VortexT,ControlPT,DragPT,NormalT] = assembly(CoordH,VortexH,ControlPH,DragPH,NormalH,CoordV,VortexV,ControlPV,DragPV,NormalV);

%Tail incidence
[CoordT,VortexT,ControlPT,DragPT,NormalT] = rotation(CoordT,VortexT,ControlPT,DragPT,NormalT,i_H,0,cr_H,x_offset_H,z_offset_H);

%Wing-body assembly
[Coord,Vortex,ControlP,DragP,Normal] = assembly(CoordW,VortexW,ControlPW,DragPW,NormalW,CoordT,VortexT,ControlPT,DragPT,NormalT);

%Plane incidence
incidence=alpha;
[Coord,Vortex,ControlP,DragP,Normal] = rotation(Coord,Vortex,ControlP,DragP,Normal,incidence,0,cr_W,x_offset_W,z_offset_W);

%Symetric plane
[Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr] = mirror (Coord,Vortex,ControlP,DragP,Normal);

%Ground-effect assembly
[Coord,Vortex,ControlP,DragP,Normal] = assembly(Coord,Vortex,ControlP,DragP,Normal,Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr);

% COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaY = [b_W/(2*Ny) b_H/(2*Ny) b_V/Ny];

Gamma = circulation(Uinf,Vortex,ControlP,Normal);
[dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala+htp+vtp');
[dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala+htp+vtp');

L = lift(dLw,dLh,dLv);
M = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala+htp+vtp');
Dind = drag(dDw,dDh,dDv);
CDparw = cdragpar(dLw,deltaY(1),Ny,cr_W,ct_W,b_W,rho,Uinf,'ala');
CDparh = cdragpar(dLh,deltaY(2),Ny,cr_H,ct_H,b_H,rho,Uinf,'htp');
CDparv = cdragpar(dLv,deltaY(3),Ny,cr_V,ct_V,b_V,rho,Uinf,'vtp');
CDpar = [CDparw CDparh CDparv];
[CL, CD, Cm] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,CDpar,M);
fprintf('Wing+VTP+HTP+Ground case: L= %f D= %f M= %f\n CL= %f CD=%f Cm=%f \n',L,Dind,M,CL,CD,Cm)

% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);
surf(Coord(:,1:2*(Ny+1),1),Coord(:,1:2*(Ny+1),2),Coord(:,1:2*(Ny+1),3));
hold on;
surf(Coord(:,2*(Ny+1)+1:4*(Ny+1),1),Coord(:,2*(Ny+1)+1:4*(Ny+1),2),Coord(:,2*(Ny+1)+1:4*(Ny+1),3));
hold on;
surf(Coord(:,4*(Ny+1)+1:5*(Ny+1),1),Coord(:,4*(Ny+1)+1:5*(Ny+1),2),Coord(:,4*(Ny+1)+1:5*(Ny+1),3));
hold on;
surf(Coord(:,5*(Ny+1)+1:7*(Ny+1),1),Coord(:,5*(Ny+1)+1:7*(Ny+1),2),Coord(:,5*(Ny+1)+1:7*(Ny+1),3));
hold on;
surf(Coord(:,7*(Ny+1)+1:9*(Ny+1),1),Coord(:,7*(Ny+1)+1:9*(Ny+1),2),Coord(:,7*(Ny+1)+1:9*(Ny+1),3));
hold on;
surf(Coord(:,9*(Ny+1)+1:10*(Ny+1),1),Coord(:,9*(Ny+1)+1:10*(Ny+1),2),Coord(:,9*(Ny+1)+1:10*(Ny+1),3));
axis equal;
