function dL = delta_lift(Vortex,Gamma,rho,Uinf)
    deltaY=(Vortex(2,:,2)-Vortex(1,:,2))';
    deltaZ=(Vortex(2,:,3)-Vortex(1,:,3))';
    delta=sqrt(deltaY.^2+deltaZ.^2);
    dL = rho*norm(Uinf)*Gamma.*delta;
end