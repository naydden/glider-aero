
% Function that computes th w_coeficcient needed to compute the induced
% Drag

function w = w_coef(Nx,Ny,vortice_mat,control,Uinf,Gamma)

N = length(vortice_mat);
w = zeros(N,1);
J = [0,1,0];

for i=1:N      % Moving along the wing drag control points
    % Coordinates of each control point
    coord_d = [control(i,1) control(i,2) control(i,3)];

    % Compute the induced velocity by all the vortex at the control
    % point (i)
    V_vortex = induced_vel_vortex(vortice_mat,coord_d);

    % Compute the total induced velocity at control point (i)
    Vind = 0;
    for l=1:N
        Vind = Vind + V_vortex(:,l) * Gamma(l);
    end
        % Compute the w coeficient
    nom = dot(cross(Uinf,Vind),J);
    denom = norm(Uinf);
    w(i) = nom/denom;
    
end

end