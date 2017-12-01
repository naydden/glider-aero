% Function that computes the induced drag of each vortex element
function [dDw,dDh,dDv] = delta_drag(Gamma,vortice_mat,control,deltaY,Nx,Ny,rho,Uinf,cas)
    dDw = zeros(Nx,2*Ny);
    dDh = zeros(Nx,2*Ny);
    dDv = zeros(Nx,2*Ny);
    w = w_coef(Nx,Ny,vortice_mat,control,Uinf,Gamma);
    switch cas
        case 'ala'
            Gamma_w = rearrange_wing(Nx,Ny,Gamma,'wing');
            w_w = rearrange_wing(Nx,Ny,w,'wing');
            for i=1:Nx
                for j=1:2*Ny
                    if i == 1
                        dDw(i,j) = rho*Gamma_w(1,j)*w_w(1,j)*deltaY(1);
                    else
                        dDw(i,j) = rho*(Gamma_w(i,j)-Gamma_w(i-1,j))*w_w(i,j)*deltaY(1);
                    end
                end
            end
        case 'ala+htp'
            Gamma_w = rearrange_wing(Nx,Ny,Gamma,'wing');
            w_w = rearrange_wing(Nx,Ny,w,'wing');
            Gamma_h = rearrange_wing(Nx,Ny,Gamma,'htp');
            w_h = rearrange_wing(Nx,Ny,w,'htp');
            for i=1:Nx
                for j=1:2*Ny
                    if i == 1
                        dDw(i,j) = rho*Gamma_w(1,j)*w_w(1,j)*deltaY(1);
                        dDh(i,j) = rho*Gamma_w(1,j)*w_h(1,j)*deltaY(2);
                    else
                        dDw(i,j) = rho*(Gamma_w(i,j)-Gamma_w(i-1,j))*w_w(i,j)*deltaY(1);
                        dDh(i,j) = rho*(Gamma_h(i,j)-Gamma_h(i-1,j))*w_h(i,j)*deltaY(2);
                    end
                end
            end            
        case 'ala+htp+vtp'
            Gamma_w = rearrange_wing(Nx,Ny,Gamma,'wing');
            w_w = rearrange_wing(Nx,Ny,w,'wing');
            Gamma_h = rearrange_wing(Nx,Ny,Gamma,'htp');
            w_h = rearrange_wing(Nx,Ny,w,'htp');
            Gamma_v = rearrange_wing(Nx,Ny,Gamma,'vtp');
            w_v = rearrange_wing(Nx,Ny,w,'vtp');
            for i=1:Nx
                for j=1:2*Ny
                    if i == 1
                        dDw(i,j) = rho*Gamma_w(1,j)*w_w(1,j)*deltaY(1);
                        dDh(i,j) = rho*Gamma_w(1,j)*w_h(1,j)*deltaY(2);
                        if(j>Ny)
                            continue
                        end
                        dDv(i,j) = rho*Gamma_v(1,j)*w_v(1,j)*deltaY(3);
                    else
                        dDw(i,j) = rho*(Gamma_w(i,j)-Gamma_w(i-1,j))*w_w(i,j)*deltaY(1);
                        dDh(i,j) = rho*(Gamma_h(i,j)-Gamma_h(i-1,j))*w_h(i,j)*deltaY(2);
                        if(j>Ny)
                            continue
                        end
                        dDv(i,j) = rho*(Gamma_v(i,j)-Gamma_v(i-1,j))*w_v(i,j)*deltaY(3);
                    end
                end
            end             
    end
end