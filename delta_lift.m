function [dLw,dLh,dLv] = delta_lift(Gamma,b,Nx,Ny,rho,Uinf,cas)
    deltaY = b/(2*Ny);
    dLw = zeros(Nx,Ny);
    dLh = zeros(Nx,Ny);
    dLv = zeros(Nx,Ny);
    switch cas
        case 'ala'
            Gamma_w = rearrange_wing(Nx,Ny,Gamma,'wing');
            for i=1:Nx
                for j=1:Ny
                    if i == 1
                        dLw(i,j) = rho*norm(Uinf)*Gamma_w(1,j)*deltaY;
                    else
                        dLw(i,j) = rho*norm(Uinf)*(Gamma_w(i,j)-Gamma_w(i-1,j))*deltaY;
                    end
                end
            end
        case 'ala+htp'
            Gamma_w = rearrange_wing(Nx,Ny,Gamma,'wing');
            Gamma_h = rearrange_wing(Nx,Ny,Gamma,'htp');
            for i=1:Nx
                for j=1:Ny
                    if i == 1
                        dLw(i,j) = rho*norm(Uinf)*Gamma_w(1,j)*deltaY;
                        dLh(i,j) = rho*norm(Uinf)*Gamma_h(1,j)*deltaY;
                    else
                        dLw(i,j) = rho*norm(Uinf)*(Gamma_w(i,j)-Gamma_w(i-1,j))*deltaY;
                        dLh(i,j) = rho*norm(Uinf)*(Gamma_h(i,j)-Gamma_h(i-1,j))*deltaY;
                    end
                end
            end            
        case 'ala+htp+vtp'
            Gamma_w = rearrange_wing(Nx,Ny,Gamma,'wing');
            Gamma_h = rearrange_wing(Nx,Ny,Gamma,'htp');
            Gamma_v = rearrange_wing(Nx,Ny,Gamma,'vtp');
            for i=1:Nx
                for j=1:Ny
                    if i == 1
                        dLw(i,j) = rho*norm(Uinf)*Gamma_w(1,j)*deltaY;
                        dLh(i,j) = rho*norm(Uinf)*Gamma_h(1,j)*deltaY;
                        dLv(i,j) = rho*norm(Uinf)*Gamma_v(1,j)*deltaY;
                    else
                        dLw(i,j) = rho*norm(Uinf)*(Gamma_w(i,j)-Gamma_w(i-1,j))*deltaY;
                        dLh(i,j) = rho*norm(Uinf)*(Gamma_h(i,j)-Gamma_h(i-1,j))*deltaY;
                        dLv(i,j) = rho*norm(Uinf)*(Gamma_v(i,j)-Gamma_v(i-1,j))*deltaY;
                    end
                end
            end             
    end
end