% Function that returns the circulation vector
function Gamma = circulation(Uinf,vortice_mat,control,n)
% vortice_mat is the matrix 5xNx3. First four rows are the A,B,C,D vertex 
% coordinates vectors of the vortex ring. Last row is 0 if vortex ring, or 1 if
% horseshoe vortex. Each column represents a single vortex.
% control is the 3xNx3
% Uinf: VECTOR of freestream velocity
N = size(vortice_mat,2); % Number of control points

a = influence_coef(vortice_mat,control,n);
% a: MATRIX of induced coefficients in the control points (NintxNint) where
% a(i,j) is theinfluence coefficient in the control point i induced by the vortex j
b = zeros(N,1);

for i = 1:N % loop over the control points
    % Calculation of the RHS coeficients (b)
    b(i) = -dot(Uinf,n(i,:));
end
% Computation of the circulation (Gamma)
Gamma = a\b;

end