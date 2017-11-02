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