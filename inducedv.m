% Function inducedv
% It calculates the induced velocity in the point P given by the vortex
% line that goes from A to B
% P = (xP, yP, zP): vector of coordinates of point P
% A = (xA, yA, zA): vector of coordinates of point A
% B = (xB, yB, zB): vector of coordinates of point B
function v = inducedv(A,B,P)
r1 = P-B;

r2 = P-A;
display(r2)
v = (abs(r1)+abs(r2))/(abs(r1)*abs(r2)*(abs(r1)*abs(r2)+dot(r1,r2)))*cross(r1,r2)/(4*pi);
end