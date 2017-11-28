
% Function that computes the induced drag of each vortex element

function dDind = delta_drag(vortice_mat,control,Gamma,b,Nx,Ny,rho,Uinf)

    deltaY = b/(2*Ny);
    dDind = zeros(Nx,Ny);
    N=Nx*Ny;
   
    w = w_coef(Nx,Ny,vortice_mat,control,Uinf,Gamma);
     
    wleft=w(1:N/2);
    wright=w(N/2+1:N);
    wMatLeft = zeros(Nx,Ny/2);
    wMatRight = zeros(Nx,Ny/2);
    
    gammaleft=Gamma(1:N/2);
    gammaright=Gamma(N/2+1:N);
    gammaMatLeft = zeros(Nx,Ny/2);
    gammaMatRight = zeros(Nx,Ny/2);
    
   
    for i = 1:Nx
        for j = 1:Ny/2
            wMatLeft(i,j)=wleft((i-1)*Ny/2+j);
            wMatRight(i,j)=wright((i-1)*Ny/2+j);
            gammaMatLeft(i,j) = gammaleft((i-1)*Ny/2+j);
            gammaMatRight(i,j) = gammaright((i-1)*Ny/2+j);
        end
    end
    
    w_new=[wMatLeft, wMatRight];
    gamma_new=[gammaMatLeft, gammaMatRight];
    
    for j=1:Ny
        for i=1:Nx
            if i == 1
                dDind(i,j) = rho*gamma_new(1,j)*w_new(1,j)*deltaY;
            else
                dDind(i,j) = rho*(gamma_new(i,j)- gamma_new(i-1,j))*w_new(i,j)*deltaY;
            end
        end
    end
    
end