% returns a Nx x Ny matrix with the inlfuence coefficients in it
function a = influence_coef(vortice_mat,control,n)
% vortice_mat is the matrix 5xNx3. First four rows are the A,B,C,D vertex 
% coordinates vectors of the vortex ring. Last row is 0 if vortex ring, or 1 if
% horseshoe vortex. Each column represents a single vortex.
% control is the  3 x N x 3 matrix with the cooordinates of each control point
% n is the matrix of normal vectors size: Nx x Ny
N = size(vortice_mat,2);
a = zeros(N,N);
for i=1:N
    % coordinates of each control point
    coord_c = [control(i,1) control(i,2) control(i,3)];
    % compute the induced velocity by all the vortex to the control
    % point (i)
    V_vortex = induced_vel_vortex(vortice_mat,coord_c);
    % now, save the influence coefficients
    for l=1:N
        n_aux = [ n(i,1) n(i,2) n(i,3)];
        a(i,l) = dot(V_vortex(:,l),n_aux);
    end
end
end