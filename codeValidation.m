% Convergence analysis to determine the number of panels to be used in the
% following calculations

% Profile geometry
m_W = 0.02;
p_W = 0.4;

% Wing geometry
cr_W = 1;
ct_W = 1;
b_W = 20;
sweep_W = 0;
dihedral_W = 0;
twist_W = 0;

% Air
alpha = 2; 
x_offset_W = 0;
z_offset_W = 0;
rho = 1.225;
Uinf = [1*cosd(alpha),0,1*sind(alpha)];

Nx = 10;
Ny = 10;
CL = 0;

deltaY = b_W/(2*Ny);
[Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
    Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);
Gamma = circulation(Uinf,Vortex,ControlP,Normal);
[dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
[dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala');
L = lift(dLw,dLh,dLv);
M = moment(dLw,dLh,dLv,Nx,Ny(i),DragP(:,:,1),'ala');
Dind = drag(dDw,dDh,dDv);
[CL, CD, Cm] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M);
fprintf('Wing validation: L= %f D= %f M= %f\n CL= %f CD=%f Cm=%f \n',L,Dind,M,CL,CD,Cm)
