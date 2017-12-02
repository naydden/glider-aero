% Function that gives the zero lift angle of attack of the wing (in the
% central section)
function alpha = ZLangle(cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset,z_offset,rho)

% Geometry
[~,Vortex,ControlP,~,Normal] = wing_assembly (cr,ct,b,...
    Nx,Ny,m,p,sweep,dihedral,twist,x_offset,z_offset);
deltaY = b/(2*Ny);

alpha0 = -twist/4; % initial value of the iteration

resta = 1;
delta = 1e-3; % error

CL = 50;
CLant = 200;

while resta>delta
    alpha = alpha0;
    CLantant = CLant;
    CLant = CL;
    
    Uinf = [1*cosd(alpha),0,1*sind(alpha)];
    
    % Computations
    Gamma = circulation(Uinf,Vortex,ControlP,Normal);
    [dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
     
    L = lift(dLw,dLh,dLv);
    [CL, ~, ~] = Coeff(cr,ct,b,Uinf,rho,L,0,0,0);
    resta = abs(CL);
    
    dalpha = 10*resta; % alpha increment
    
    % new alpha
    if(CL>0)
        alpha0 = alpha0-dalpha;
    else
        alpha0 = alpha0+dalpha;
    end
    
    % In order to avoid an infinite loop
    if CL==CLantant
        resta = 0;
    end
    
end
