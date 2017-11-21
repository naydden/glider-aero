% returns a matrix with the induced velocity at each control point
% size of matrix: Nx x Ny x 3
function V = induced_vel(vort_p,con_p)
% vort_p is the 3x(Nx+1)x(Ny+1) matrix
% con_p is the 3x Nx x Ny
Nx = size(con_p);
Nx = Nx(1);
Ny = Nx(2);
k=1;
for i=1:Nx
    for j=1:Ny
        coord_p = [con_p(i,j,1) con_p(i,j,2) con_p(i,j,3)];
        V_vortx = induced_vel_vortex(vort_p,coord_p)
% Pensar
%         for l=1:Nx*Ny
%             a(i,j) = V_vortx(:,l)*n(i,j)
%         end
    end
end
end