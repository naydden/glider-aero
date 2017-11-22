
% Function that computes the induced drag of each vortex element

function dDind = delta_drag(vortice_mat,control,Gamma,b,Nx,Ny,ro,Uinf)

    deltaY = b/(2*Ny);
    dDind = zeros(Nx,Ny);
    Gamma_new = zeros(Nx,Ny);
    
    w = w_coef(vortice_mat,control,Uinf,Gamma)
    
    for j=1:Ny
        for i=1:Nx
            
            % Change from vector to matrix
            Gamma_new(i,j) = Gamma((i-1)*Ny+j);
            
            if i == 1
                dDind(i,j) = ro*Gamma_new(1,j)*w(1,j)*deltaY;
                
            else
                dDind(i,j) = ro*(Gamma_new(i,j)- Gamma_new(i-1,j))*w(i,j)*deltaY;
            end
            
        end
    end
    
end