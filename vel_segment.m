% function that returns induced velocity by a line vortex b-c, at any point
% horseshoe vortex:
% B-----------C
% |           |
% |           |
% |           |
% |           |
% A           D
function V = vel_segment (x_b,y_b,z_b,x_c,y_c,z_c,x,y,z)
    r1 = [x-x_b; y-y_b; z-z_b];
    r2 = [x-x_c; y-y_c; z-z_c];
    V = 1/(4*pi)*(norm(r1)+norm(r2))/(norm(r1)*norm(r2)*(norm(r1)*norm(r2) + r1'*r2))*cross(r1,r2);
end