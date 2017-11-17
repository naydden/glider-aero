%function that returns induced velocity by a given vortex line a-b
% at point x,y,z
function V = vel_segment (x_a,y_a,z_a,x_b,y_b,z_b,x,y,z,cas)
switch cas
    case 'line'
        r1 = [x-x_a; y-y_a; z-z_a]
        r2 = [x-x_b; y-y_b; z-z_b]
        V = 1/(4*pi)*(norm(r1)+norm(r2))/(norm(r1)*norm(r2)*(norm(r1)*norm(r2) + r1'*r2))*cross(r1,r2)
    case 'semi-line'
        ur = [-1 0 0]
        r2 = [x-x_b; y-y_b; z-z_b]
        ur2 = r2/norm(r2);
        V = 1/(4*pi)*cross(ur,r2)/(norm(ur,r2)^2)*(1-ur'*ur2);   
    case 'semi-line-e'

end