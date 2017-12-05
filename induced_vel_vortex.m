% returns the induced velocity to a point by all the vortex
function V = induced_vel_vortex(vort_p,con_p)
% vort_p is the (Nx+1)x(Ny+1)x3 matrix
% con_p is x y z of the control point of the vortex we are analyzing
N = size(vort_p,2);

x = con_p(1);
y = con_p(2);
z = con_p(3);

V = [0; 0; 0;];
for i=1:N
        x_b = vort_p(1,i,1); y_b = vort_p(1,i,2); z_b = vort_p(1,i,3);
        x_c = vort_p(2,i,1); y_c = vort_p(2,i,2); z_c = vort_p(2,i,3);
        V_ab = vel_segment(x_b,y_b,z_b,x_c,y_c,z_c,x,y,z);
        AB=[x_c-x_b,y_c-y_b,z_c-z_b];
        AP=[x-x_b,y-y_b,z-z_b];
        if abs( abs(dot(AB,AP)) - norm(AB)*norm(AP) ) < 1e-7
            V_ab = [0; 0; 0;];
        end
%         if isnan(V_ab(1)) || isinf(V_ab(1)) || isnan(V_ab(2)) || isinf(V_ab(2)) || isnan(V_ab(3)) || isinf(V_ab(3))
%             V_ab = [0; 0; 0;];
%         end
        % 0 equals vortex ring
        if vort_p(5,i,1) == 0          
            x_b = vort_p(2,i,1); y_b = vort_p(2,i,2); z_b = vort_p(2,i,3);
            x_c = vort_p(3,i,1); y_c = vort_p(3,i,2); z_c = vort_p(3,i,3);
            V_bc = vel_segment(x_b,y_b,z_b,x_c,y_c,z_c,x,y,z);
            x_b = vort_p(3,i,1); y_b = vort_p(3,i,2); z_b = vort_p(3,i,3);
            x_c = vort_p(4,i,1); y_c = vort_p(4,i,2);  z_c = vort_p(4,i,3);
            V_cd = vel_segment(x_b,y_b,z_b,x_c,y_c,z_c,x,y,z);
            CD=[x_c-x_b,y_c-y_b,z_c-z_b];
            CP=[x-x_b,y-y_b,z-z_b];
            if abs( abs(dot(CD,CP)) - norm(CD)*norm(CP) ) < 1e-7
                V_cd = [0; 0; 0;];
            end    
%             if isnan(V_cd(1)) || isinf(V_cd(1)) || isnan(V_cd(2)) || isinf(V_cd(2)) || isnan(V_cd(3)) || isinf(V_cd(3))
%                 V_cd = [0; 0; 0;];
%             end
            x_b = vort_p(4,i,1); y_b = vort_p(4,i,2); z_b = vort_p(4,i,3);
            x_c = vort_p(1,i,1); y_c = vort_p(1,i,2); z_c = vort_p(1,i,3);
            V_da = vel_segment(x_b,y_b,z_b,x_c,y_c,z_c,x,y,z);              
        else
            x_e = vort_p(2,i,1);
            y_e = vort_p(2,i,2);
            z_e = vort_p(2,i,3);
            V_bc = vel_seminfline(x_e,y_e,z_e,x,y,z,'CD');
            V_cd = [0; 0; 0;];
            x_e = vort_p(1,i,1);
            y_e = vort_p(1,i,2);
            z_e = vort_p(1,i,3);
            V_da = vel_seminfline(x_e,y_e,z_e,x,y,z,'AB');            
        end
        V(:,i) = V_ab + V_bc + V_cd + V_da;
end
end