
%% Function that returns the lift coeficient 

function CL = lift_coeff(dLw,dLh,Sw_S,Sh_S,rho,Uinf,Nx,Ny)

    Lw = sum(sum(dLw));
    Lh = sum(sum(dLh));
    
    CLw = 2*Lw / ( rho*Sw_S*norm(Uinf)^2 );
    CLh = 2*Lh / ( rho*Sh_S*norm(Uinf)^2 );
    
    CL = CLw+CLh;

end