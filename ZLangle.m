% Function that gives the zero lift angle of attack of the wing (in the
% central section)
function alpha = ZLangle(cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset,z_offset,rho)

alpha0 = 0;
resta = 1;
delta = 1e-3; % error

while resta>delta
    alpha = alpha0;
    Uinf= [1*cosd(alpha),0,1*sind(alpha)];
    [CL, ~] = Coeff(cr,ct,b,Nx,Ny,m,p,sweep,dihedral,twist,x_offset,z_offset,Uinf,rho);
    resta = abs(CL);
    
    % new alpha
    if(CL>0)
        alpha0 = alpha0-0.01;
    else
        alpha0 = alpha0+0.01;
    end
    
end
display(CL)
    