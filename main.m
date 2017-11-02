%% Solving the aerodynamics of a conventional glider
% Authors:
%   - Silvia Gonzalez
%   - Laura Pla
%   - Josep Maria Serra Moncunill
%   - Boyan Naydenov
% Subject: Aerdonyamics, Flight and Orbital Mechanics.
% Date: December 5th, 2017
%% Data
% St = horizontal tail (HTP) area
% Sv = vertical tail (VTP) area
% S = wing area
% lt = distance given in drawing
% Cm = mean aerodynamic chord of the wing
St_S = 1/8;
Sv_S = 2/3;
lt_Cm = 4;
% Assumptions:
% AR of wing and HTP is high enough (>6) -> lifting line can be assumed
% Ct/Cr (tip to chord) ratio is left free.
%% Wing:NACA 2412 airfoil
i_w0 = 0; %the chord of the center section of the wing has zero incidence 
% angle with respect to the rod.
%% HTP NACA 0009
% symmetric airfoil
i_w_htp_0 = -4; %degrees
%% VTP NACA 0009
% symmetric airfoil

%% Possible methods
% - Vortex lattice
% - Lifting line

%% ONLY FOR ISOLATED WING
%% Part 1: Compute ZL angle of wing for twist (0 to 8 deg) and CD.
%% Part 2: PLotting wing's aerodynamic polar for alpha 0 to 10 deg.
%% Part 3: Assumption -> Ground Effect. Plot CL and CD for alpha 6deg, against AR 0.075Ao to 1.25Ao
% Ao is the nominal aspect ratio


%% WING + HTP + VTP
%% Part 4: CL and CD and Xcm for Sum(M)=0 when alpha = 6deg
%% Part 5: Assumption -> ground effect. CL, CD and CM_cm