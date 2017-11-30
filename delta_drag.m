
% Function that computes the induced drag of each vortex element

function dDind = delta_drag(vortice_mat,control,Gamma,deltaY,Nx,Ny,rho,Uinf)

    w = w_coef(Nx,Ny,vortice_mat,control,Uinf,Gamma);
  
    n = length(Gamma);
    Nw = ceil(n/(Nx*Ny));
    
    if Nw>5 | Nw==4
        n=n/2; Nw=Nw/2;
    end
    
    w_new = zeros(Nx,Ny*Nw);
    Gamma_new = zeros(Nx,Ny*Nw);
    dDind = zeros(Nx,Ny*Nw);  
    
    for Nw = 1:ceil(n/(2*Nx*Ny));
        for i = 1:Nx
            if Nw == 3  %Cas de l'estabilitzador vertical
                for j = 1:Ny
                    w_new(i,j+2*Ny*(Nw-1)) = w((i-1)*Ny + j + 4*Ny*(Nw-1));
                    Gamma_new(i,j+2*Ny*(Nw-1)) = Gamma((i-1)*Ny + j + 4*Ny*(Nw-1));
                    
                    if i == 1
                        dDind(i,j+2*Ny*(Nw-1)) = rho*Gamma_new(i,j+2*Ny*(Nw-1))*w_new(1,j)*deltaY(Nw);
                    else
                        dDind(i,j+2*Ny*(Nw-1)) = rho*(Gamma_new(i,j+2*Ny*(Nw-1))- Gamma_new(i-1,j+2*Ny*(Nw-1)))*w_new(i,j+2*Ny*(Nw-1))*deltaY(Nw);
                    end
                    
                end
            else   % Cas ala i estabilitzador horitzontal
                for j = 1:2*Ny
                    w_new(i,j+2*Ny*(Nw-1)) = w((i-1)*2*Ny + j + 4*Ny*(Nw-1));
                    Gamma_new(i,j+2*Ny*(Nw-1)) = Gamma((i-1)*2*Ny + j + 4*Ny*(Nw-1)); 
                   
                    if i == 1
                        dDind(i,j+2*Ny*(Nw-1)) = rho*Gamma_new(i,j+2*Ny*(Nw-1))*w_new(1,j)*deltaY(Nw);
                    else
                        dDind(i,j+2*Ny*(Nw-1)) = rho*(Gamma_new(i,j+2*Ny*(Nw-1))- Gamma_new(i-1,j+2*Ny*(Nw-1)))*w_new(i,j+2*Ny*(Nw-1))*deltaY(Nw);
                    end
                    
                end
            end
        end
    end
    
end