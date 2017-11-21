% returns the induced velocity to a point by all the vortex
function V = induced_vel_vortex(vort_p,con_p)
% vort_p is the (Nx+1)x(Ny+1)x3 matrix
% con_p is x y z of the control point of the vortex we are analyzing
Nx = size(vort_p,1)-1;
Ny = size(vort_p,2)-1;

x = con_p(1);
y = con_p(2);
z = con_p(3);

V = [0; 0; 0;];
k = 1;
for i=1:Nx
    for j=1:Ny
        x_b = vort_p(i,j,1);
        y_b = vort_p(i,j,2);
        z_b = vort_p(i,j,3);
        x_c = vort_p(i,j+1,1);
        y_c = vort_p(i,j+1,2);
        z_c = vort_p(i,j+1,3);
        V_ab = vel_segment(x_b,y_b,z_b,x_c,y_c,z_c,x,y,z);
        if i~=Nx            
            x_b = vort_p(i,j+1,1);
            y_b = vort_p(i,j+1,2);
            z_b = vort_p(i,j+1,3);
            x_c = vort_p(i+1,j+1,1);
            y_c = vort_p(i+1,j+1,2);
            z_c = vort_p(i+1,j+1,3);
            V_bc = vel_segment(x_b,y_b,z_b,x_c,y_c,z_c,x,y,z);

            x_b = vort_p(i+1,j+1,1);
            y_b = vort_p(i+1,j+1,2);
            z_b = vort_p(i+1,j+1,3);
            x_c = vort_p(i+1,j,1);
            y_c = vort_p(i+1,j,2);
            z_c = vort_p(i+1,j,3);
            V_cd = vel_segment(x_b,y_b,z_b,x_c,y_c,z_c,x,y,z);

            x_b = vort_p(i+1,j,1);
            y_b = vort_p(i+1,j,2);
            z_b = vort_p(i+1,j,3);
            x_c = vort_p(i,j,1);
            y_c = vort_p(i,j,2);
            z_c = vort_p(i,j,3);
            V_da = vel_segment(x_b,y_b,z_b,x_c,y_c,z_c,x,y,z);              
        else
            x_e = vort_p(i,j+1,1);
            y_e = vort_p(i,j+1,2);
            z_e = vort_p(i,j+1,3);
            V_bc = vel_seminfline(x_e,y_e,z_e,x,y,z,'CD');
            V_cd = [0; 0; 0;];
            x_e = vort_p(i,j,1);
            y_e = vort_p(i,j,2);
            z_e = vort_p(i,j,3);
            V_da = vel_seminfline(x_e,y_e,z_e,x,y,z,'AB');            
        end
        V(:,k) = V_ab + V_bc + V_cd + V_da;
        k = k + 1;
    end
end
end