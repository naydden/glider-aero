
%% Solving the aerodynamics of a conventional glider using the vortex lattice method
% Authors:
%   - Silvia Gonzalez
%   - Laura Pla
%   - Josep Maria Serra Moncunill
%   - Boyan Naydenov
% Subject: Aerodynamics, Flight and Orbital Mechanics.
% Date: December 5th, 2017

%%
clear
clc


%% Data

% St = horizontal tail (HTP) area
% Sv = vertical tail (VTP) area
% S = wing area
% lt = distance given in drawing
% cm = mean aerodynamic chord of the wing
% m = maximum chamber
% p = location of maximum chamber (% of the chord)
% i = incidence angle

% Wing:NACA 2412 airfoil
m_W = 0.02;
p_W = 0.4;
i_W = 0; %the chord of the center section of the wing has zero incidence angle with respect to the rod

% HTP NACA 0009
m_H = 0;
p_H = 0;
i_H = 0; % degrees

% VTP NACA 0009
m_V = 0;
p_V = 0;

St_S = 1/8;
Sv_S = 2/3;
lt_cm = 4;

rho = 1.225; % density [kg/m^3]

Nx=5; Ny=10;

% Assumptions:
% AR of wing and HTP is high enough (>6) -> lifting line can be assumed
% Ct/Cr (tip to chord) ratio is left free.

%% Part 1: Compute ZL angle of wing for twist (0 to 8 deg) and CD.

% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0; % angle of attack [degrees]
Uinf = [1*cosd(alpha),0,1*sind(alpha)];

% Wing
cr_W=1; ct_W=1*cr_W; b_W=10;
sweep_W=0; dihedral_W=0; twist_W=0; 
x_offset_W=0; z_offset_W=0;

[Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
    Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);

% COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaY = b_W/(2*Ny);

Gamma = circulation(Uinf,Vortex,ControlP,Normal);
[dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
[dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala');

L = lift(dLw,dLh,dLv);
M = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala');
Dind = drag(dDw,dDh,dDv);
[CL, CD, Cm] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M);
fprintf('Wing case: L= %f D= %f M= %f\n CL= %f CD=%f Cm=%f \n',L,Dind,M,CL,CD,Cm)

% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
surf(Coord(:,1:2*(Ny+1),1),Coord(:,1:2*(Ny+1),2),Coord(:,1:2*(Ny+1),3));
axis equal;

%% Part 2: PLotting wing's aerodynamic polar for alpha 0 to 10 deg.


%% Part 3: Assumption -> Ground Effect. 
% Plot CL and CD for alpha 6deg, against AR 0.075Ao to 1.25Ao
% Ao is the nominal aspect ratio

% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 3;
Uinf = [1,0,0];

% Wing
cr_W=1; ct_W=1*cr_W; b_W=10;
sweep_W=0; dihedral_W=0; twist_W=0; 
x_offset_W=0; z_offset_W=1;

[Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,Nx,...
    Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);

%Plane incidence
incidence=alpha;
[Coord,Vortex,ControlP,DragP,Normal] = rotation(Coord,Vortex,ControlP,DragP,Normal,incidence,0,cr_W,x_offset_W,z_offset_W);

%Symetric plane
[Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr] = mirror (Coord,Vortex,ControlP,DragP,Normal);

%Ground-effect assembly
[Coord,Vortex,ControlP,DragP,Normal] = assembly(Coord,Vortex,ControlP,DragP,Normal,Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr);

% COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaY = b_W/(2*Ny);

Gamma = circulation(Uinf,Vortex,ControlP,Normal);
[dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
[dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala');

L = lift(dLw,dLh,dLv);
M = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala');
Dind = drag(dDw,dDh,dDv); 
[CL, CD, Cm] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M);
fprintf('Wing case + ground: L= %f D= %f M= %f\n CL= %f CD=%f Cm=%f \n',L,Dind,M,CL,CD,Cm)

% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);
surf(Coord(:,1:2*(Ny+1),1),Coord(:,1:2*(Ny+1),2),Coord(:,1:2*(Ny+1),3));
hold on;
surf(Coord(:,2*(Ny+1)+1:4*(Ny+1),1),Coord(:,2*(Ny+1)+1:4*(Ny+1),2),Coord(:,2*(Ny+1)+1:4*(Ny+1),3));
axis equal;


%% Part 4: CL and CD and Xcm for Sum(M)=0 when alpha = 6deg

% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0;
Uinf = [1*cosd(alpha),0,1*sind(alpha)];

% Wing
cr_W=1; ct_W=1*cr_W; b_W=10;
sweep_W=0; dihedral_W=0; twist_W=0; 
x_offset_W=0; z_offset_W=0;

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

% COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaY = [b_W/(2*Ny) b_H/(2*Ny) b_V/Ny];

Gamma = circulation(Uinf,Vortex,ControlP,Normal);
[dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala+htp+vtp');
[dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala+htp+vtp');

L = lift(dLw,dLh,dLv);
M = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala+htp+vtp');
Dind = drag(dDw,dDh,dDv);
[CL, CD, Cm] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M);
fprintf('Wing+VTP+HTP case: L= %f D= %f M= %f\n CL= %f CD=%f Cm=%f \n',L,Dind,M,CL,CD,Cm)

% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);
surf(Coord(:,1:2*(Ny+1),1),Coord(:,1:2*(Ny+1),2),Coord(:,1:2*(Ny+1),3));
hold on;
surf(Coord(:,2*(Ny+1)+1:4*(Ny+1),1),Coord(:,2*(Ny+1)+1:4*(Ny+1),2),Coord(:,2*(Ny+1)+1:4*(Ny+1),3));
hold on;
surf(Coord(:,4*(Ny+1)+1:5*(Ny+1),1),Coord(:,4*(Ny+1)+1:5*(Ny+1),2),Coord(:,4*(Ny+1)+1:5*(Ny+1),3));
axis equal;

%% Part 5: Assumption -> ground effect. CL, CD and CM_cm

% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 3;
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
[CL, CD, Cm] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M);
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
