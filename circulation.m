% Function that returns the circulation vector
function Gamma = circulation(Uinf,vortice_mat,control,n)
% a: MATRIX of induced coefficients in the control points (NintxNint) where
% a(i,j) is theinfluence coefficient in the control point i induced by the vortex j
% Uinf: VECTOR of freestream velocity
Nx = size(vortice_mat,1)-1;
Ny = size(vortice_mat,2)-1;
Nint = Nx*Ny; % Number of control points

a = influence_coef(vortice_mat,control,n);
b = zeros(Nint);
n_vector = zeros(Nx*Ny,3);
% transformation from Matrix to vector of n
for l = 1:Nx
    for m = 1:Ny
        n_vector((l-1)*Ny+m,:) = n(l,m,:);
    end
end

for i = 1:Nint % loop over the control points
    % Calculation of the RHS coeficients (b)
    b(i) = -dot(Uinf,n_vector(i,:));
end

% Computation of the circulation (Gamma)
Gamma = a\b;

end