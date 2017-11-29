% Convergence analysis to determine the number of panels to be used in the
% following calculations

% Profile geometry
m_w = 0.02;
p_w = 0.4;

% Wing geometry
cr = 1;
ct = 1;
b = 20;
sweep = 0;
dihedral = 0;
twist = 0;

% Air
alpha = 2; 
x_offset = 0;
z_offset = 0;
rho = 1.225;
Uinf = [1*cosd(alpha),0,1*sind(alpha)];

%% Ny

Nx = 10;
Ny = 10:10:20;
CL = zeros(1,size(Ny,2));

for i = 1:size(Ny,2)
    [CL(i), ~] = Coeff(cr,ct,b,Nx,Ny(i),m_w,p_w,sweep,dihedral,twist,x_offset,z_offset,Uinf,rho);
    display(Ny(i));
end

figure(1);
plot(Ny,CL);
xlabel('N_{y}');
ylabel('c_{L}');
title('N_{x} = 10');
grid on;

%% Nx

Numberx = 2:2:10;
Numbery = 3;
CL = zeros(1,size(Numberx,2));

for i = 1:size(Numberx,2)
    [CL(i), ~] = Coeff(cr,ct,b,Numberx(i),Numbery,m_w,p_w,sweep,dihedral,twist,x_offset,z_offset,Uinf,rho);
    display(Numberx(i));
end

figure(20);  
plot(Numberx,CL);
xlabel('N_{x}');
ylabel('c_{L}');
title('N_{y} = 40');
grid on;