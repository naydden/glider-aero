% Convergence analysis to determine the number of panels to be used in the
% following calculations

% Profile geometry
m_W = 0.02;
p_W = 0.4;

% Wing geometry
lambda = 0.3;
A_ratio = 26;
cr_W = 1;
ct_W = lambda*cr_W;
b_W=A_ratio*0.5*(cr_W+ct_W);
sweep_W = 0;
dihedral_W = 0;
twist_W = 0;

% Air
alpha = 2; 
x_offset_W = 0;
z_offset_W = 0;
rho = 1.225;
Uinf = [1*cosd(alpha),0,1*sind(alpha)];
CDpar = [0 0 0];

%% Ny

Nx = 10;
Ny = 10:5:50;
CL = zeros(1,size(Ny,2));

for i = 1:size(Ny,2)
    deltaY = b_W/(2*Ny(i));
    [Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
        Nx,Ny(i),m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);
    Gamma = circulation(Uinf,Vortex,ControlP,Normal);
    [dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny(i),rho,Uinf,'ala');
    L = lift(dLw,dLh,dLv);
    [CL(i), ~, ~] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,0,CDpar,0);
    display(Ny(i));
end

figure(1);
plot(Ny,CL);
xlabel('N_{y}');
ylabel('c_{L}');
title(['N_{x} = ' num2str(Nx)]);
grid on;

%% Nx

Numberx = 2:2:12;
Numbery = 30;
CLift = zeros(1,size(Numberx,2));
deltaY = b_W/(2*Numbery);

for i = 1:size(Numberx,2)
    [Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
        Numberx(i),Numbery,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);
    Gamma = circulation(Uinf,Vortex,ControlP,Normal);
    [dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Numberx(i),Numbery,rho,Uinf,'ala');
    L = lift(dLw,dLh,dLv);
    [CLift(i), ~, ~] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,0,CDpar,0);
    display(Numberx(i));
end

figure(2);  
plot(Numberx,CLift);
xlabel('N_{x}');
ylabel('c_{L}');
title(['N_{y} = ' num2str(Numbery)]);
grid on;