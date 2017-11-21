function dL = delta_lift(Gamma,b,Nx,Ny,ro,Uinf)
    deltaY = (2*Ny)/b;
    dL = zeros(Nx,Ny);
    Gamma_new = zeros(Nx,Ny);
    for j=1:Ny
        for i=1:Nx
            % change from vector to matrix
            Gamma_new(i,j) = Gamma((i-1)*Ny+j);
            % Not sure if what Uinf to use. Using norm for now! ALERT!
            if i == 1
                dL(i,j) = ro*norm(Uinf)*Gamma_new(1,j)*deltaY;
            else
                dL(i,j) = ro*norm(Uinf)*(Gamma_new(i,j)-Gamma_new(i-1,j))*deltaY;
            end
        end
    end
end