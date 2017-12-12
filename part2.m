
% Wing
twist_W=0;

alpha = 0:1:10;
miau = size(alpha,2);
CL = zeros(1,miau);
CD = zeros(1,miau);

% Geometry
[Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
    Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);
deltaY = b_W/(2*Ny);

for i = 1:miau
    
    Uinf = [1*cosd(alpha(i)),0,1*sind(alpha(i))];
    Gamma = circulation(Uinf,Vortex,ControlP,Normal);
    [dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
    [dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala');
    
    L = lift(dLw,dLh,dLv);
    Dind = drag(dDw,dDh,dDv);
    CDparw = cdragpar(dLw,deltaY,Ny,cr_W,ct_W,b_W,rho,Uinf,'ala');
    CDpar = [CDparw 0 0];
    [CL(i), CD(i), ~] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,CDpar,0);

end

figure(3);
plot(CL,CD);
xlabel('C_{L}')
ylabel('C_{D}')
grid on;

figure(5);
surf(Coord(:,1:2*(Ny+1),1),Coord(:,1:2*(Ny+1),2),Coord(:,1:2*(Ny+1),3));
axis equal;