function dL = delta_lift(Gamma,b,Nx,Ny,rho,Uinf)
    deltaY = b/(2*Ny);
    dLift = rho*norm(Uinf)*Gamma*deltaY;
    N=Nx*Ny;
    
    dLleft=dLift(1:N/2);
    dLright=dLift(N/2+1:N);
    
    dLMatLeft = zeros(Nx,Ny/2);
    dLMatRight = zeros(Nx,Ny/2);
    
    for i = 1:Nx
        for j = 1:Ny/2
            dLMatLeft(i,j)=dLleft((i-1)*Ny/2+j);
            dLMatRight(i,j)=dLright((i-1)*Ny/2+j);
        end
    end
    dL=[dLMatLeft, dLMatRight];
end