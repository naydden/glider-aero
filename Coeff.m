% Function that calculates the aerodynamic coefficients of the wing (CLift
% and CDrag)
function [CL, CD] = Coeff(cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset,z_offset,Uinf,rho)

% Geometry
[~,Vortex,ControlP,DragP,Normal] = wing_assembly (cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset,z_offset);

Gamma = circulation(Uinf,Vortex,ControlP,Normal);

% Lift
dL = delta_lift(Vortex,Gamma,rho,Uinf);
L = lift(dL,0);

% Drag
dDind = delta_drag(Vortex,DragP,Gamma,b,Nx,Ny,rho,Uinf);
Dind = drag(dDind,Nx,Ny);

% Coefficients
S = b*(cr+ct)/2;
CL = 2*L/(rho*norm(Uinf)^2*S);
CD = 2*Dind/(rho*norm(Uinf)^2*S);

end