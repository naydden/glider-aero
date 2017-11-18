% function that returns induced velocity by a semi-infinite line vortex, at any point
% horseshoe vortex
% B-----------C
% |           |
% |           |
% |           |
% |           |
% A           D         
function V = vel_seminfline(x_e,y_e,z_e,x,y,z,cas)
        if(cas == 'AB')
            % case AB semiline
            ur = [-1; 0; 0];
            r2 = [x-x_e; y-y_e; z-z_e];
            ur2 = r2/norm(r2);
            V = 1/(4*pi)*cross(ur,r2)/(norm(cross(ur,r2))^2)*(1-ur'*ur2);
        elseif(cas == 'CD')
            % case CD semiline
            ur = [1; 0; 0];
            r1 = [x-x_e; y-y_e; z-z_e];
            ur1 = r1/norm(r1);
            V = 1/(4*pi)*cross(ur,r1)/(norm(cross(ur,r1))^2)*(ur'*ur1+1);
        end
end