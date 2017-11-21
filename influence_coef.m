% returns a Nx x Ny matrix with the inlfuence coefficients in it
function a = influence_coef(vortice_mat,control,n)
% vortice_mat is the (Nx+1)x(Ny+1)x3 matrix with the coordinates of each vortex ring 
% control is the  Nx x Ny x 3 matrix with the cooordinates of each control point
% n is the matrix of normal vectors size: Nx x Ny
Nx = size(vortice_mat,1)-1;
Ny = size(vortice_mat,2)-1;
a = zeros(Nx*Ny,Nx*Ny);
k=1;
for i=1:Nx
    for j=1:Ny
        % coordinates of each control point
        coord_c = [control(i,j,1) control(i,j,2) control(i,j,3)];
        % compute the induced velocity by all the vortex to the control
        % point (i,j)
        V_vortex = induced_vel_vortex(vortice_mat,coord_c);
        % now, save the influence coefficients
        for l=1:Nx*Ny
            n_aux = [ n(i,j,1) n(i,j,2) n(i,j,3)];
            a(k,l) = dot(V_vortex(:,l),n_aux);
        end
        k = k + 1;
    end
end
end