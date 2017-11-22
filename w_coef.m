
% Function that computes th w_coeficcient needed to compute the induced
% Drag

function w = w_coef(vortice_mat,control,Uinf,Gamma)

Nx = size(vortice_mat,1)-1;
Ny = size(vortice_mat,2)-1;
w = zeros(Nx,Ny);
J = [0,1,0];

for i=1:Nx
    for j=1:Ny
        
        % Coordinates of each control point
        coord_c = [control(i,j,1) control(i,j,2) control(i,j,3)];
        
        % Compute the induced velocity by all the vortex at the control
        % point (i,j)
        V_vortex = induced_vel_vortex(vortice_mat,coord_c);
        
        % Compute the total induced velocity at control point (i,j)
        for l=1:Nx*Ny
            Vind = Vind + V_vortex(:,l) * Gamma(l);
        end
        
        % Compute the w coeficient
        nom = dot(cross(Uinf,Vind),J);
        denom = norm(Uinf);
        w(i,j) = nom/denom;
        
    end
end

end