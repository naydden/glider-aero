% Function that gives the zero lift angle of attack of the wing (in the
% central section)
function alpha = ZLangle(cr_W,ct_W,b_W,Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist,x_offset_W,z_offset_W,rho)

% Geometry
[~,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
    Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);

alpha0 = -twist/4; % initial value of the iteration
dalpha = 0.01; % alpha increment

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
    [dLw,dLh,dLv] = delta_lift(Gamma,b_W,Nx,Ny,rho,Uinf,'ala');
    dDind = delta_drag(Vortex,DragP,Gamma,deltaY,Nx,Ny,rho,Uinf);
    
    L = lift(dLw,dLh,dLv);
    M = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala');
    Dind = drag(dDind,Nx,Ny);
    [CL, ~, ~] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M);
    resta = abs(CL);
    
    
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
display(CL)
    