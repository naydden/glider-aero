% Function that returns the circulation vector
function Gamma = circulation(v, n, Uinf)
% v: MATRIX of velocities induced in the control points (NintxNint) where
% v(i,j) is the velocity in the control point i induced by the vortex j
% n: MATRIX of vectors normal to the ZLL in the control point (3xNint)
% Uinf: VECTOR of freestream velocity

Nint = size(n); % Number of control points

a = zeros(Nint,Nint);
b = zeros(Nint);

for i = 1:Nint % loop over the control points
    for j = 1:Nint % loop over the vortex points
        % Calculation of the influence coefficients (a)
        a(i,j) = dot(v(i,j),n(:,i));
    end
    % Calculation of the RHS coeficients (b)
    b(i) = -dot(Uinf,n(:,i));
end

% Computation of the circulation (Gamma)
Gamma = a\b;

end