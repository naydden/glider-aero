%% Solving the aerodynamics of a conventional glider
% Authors:
%   - Silvia Gonzalez
%   - Laura Pla
%   - Josep Maria Serra Moncunill
%   - Boyan Naydenov
% Subject: Aerodynamics, Flight and Orbital Mechanics.
% Date: December 5th, 2017
%% Data
% St = horizontal tail (HTP) area
% Sv = vertical tail (VTP) area
% S = wing area
% lt = distance given in drawing
% cm = mean aerodynamic chord of the wing
St_S = 1/8;
Sv_S = 2/3;
lt_cm = 4;
% Assumptions:
% AR of wing and HTP is high enough (>6) -> lifting line can be assumed
% Ct/Cr (tip to chord) ratio is left free.
%% Wing:NACA 2412 airfoil
m_w = 0.02;
p_w = 0.4;
i_w0 = 0; %the chord of the center section of the wing has zero incidence 
% angle with respect to the rod.
%% HTP NACA 0009
% symmetric airfoil
m_htp = 0;
p_htp = 0;
i_w_htp_0 = -4; %degrees
%% VTP NACA 0009
% symmetric airfoil

%% Possible methods
% - Vortex lattice
% - Lifting line
%% ONLY FOR ISOLATED WING
%% Input
cr=1; ct=1; b=20; Nx=3; Ny=3; sweep=0; dihedral=0; twist=0; alpha=0; 
ro = 1.225; Uinf= [1*cosd(alpha),0,1*sind(alpha)];
%% Preliminary
[Coord,CoordP,CoordC,CoordD,n] = geometry (cr,ct,b,Nx,Ny,m_w,p_w,sweep,dihedral,twist);
Gamma = circulation(Uinf,CoordP,CoordC,n);
dL = delta_lift(Gamma,b,Nx,Ny,ro,Uinf);
dDind = delta_drag(CoordP,CoordD,Gamma,b,Nx,Ny,ro,Uinf);
L = lift(dL,Nx,Ny);
M = moment(dL,Nx,Ny,CoordP(:,:,1));
Dind = drag(dDind,Nx,Ny);
%% Part 1: Compute ZL angle of wing for twist (0 to 8 deg) and CD.
%% Part 2: PLotting wing's aerodynamic polar for alpha 0 to 10 deg.
%% Part 3: Assumption -> Ground Effect. Plot CL and CD for alpha 6deg, against AR 0.075Ao to 1.25Ao
% Ao is the nominal aspect ratio


%% WING + HTP + VTP
%% Part 4: CL and CD and Xcm for Sum(M)=0 when alpha = 6deg
%% Part 5: Assumption -> ground effect. CL, CD and CM_cm