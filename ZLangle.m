% Function that gives the zero lift angle of attack of the wing (in the
% central section)
function alpha = ZLangle(cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset,z_offset,rho)

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
    
    Uinf= [1*cosd(alpha),0,1*sind(alpha)];
    [CL, ~] = Coeff(cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset,z_offset,Uinf,rho);
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
    