% Function that calculates the aerodynamic coefficients of the wing (CLift
% and CDrag)
function [CL, CD, Cm] = Coeff(cr,ct,b,Uinf,rho,L,Dind,M)

% Geometry
S = b*(cr+ct)/2;
MGC = S/b;

q = 0.5*rho*norm(Uinf)^2;

% Coefficients
CL = L/(q*S);
% CD = cdrag('2412',CL)+Dind/(q*S);
CD = Dind/(q*S);
Cm = M/(q*S*MGC);

end