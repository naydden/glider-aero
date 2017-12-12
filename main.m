
%% Solving the aerodynamics of a conventional glider using the vortex lattice method
% Authors:
%   - Silvia Gonzalez
%   - Laura Pla
%   - Josep Maria Serra Moncunill
%   - Boyan Naydenov
% Subject: Aerodynamics, Flight and Orbital Mechanics.
% Date: December 5th, 2017

%%
clear
clc

%% Data

% St = horizontal tail (HTP) area
% Sv = vertical tail (VTP) area
% S = wing area
% lt = distance given in drawing
% cm = mean aerodynamic chord of the wing
% m = maximum chamber
% p = location of maximum chamber (% of the chord)
% i = incidence angle

% Wing:NACA 2412 airfoil
m_W = 0.02;
p_W = 0.4;
i_W = 0; %the chord of the center section of the wing has zero incidence angle with respect to the rod

% HTP NACA 0009
m_H = 0;
p_H = 0;
i_H = 0; % degrees

% VTP NACA 0009
m_V = 0;
p_V = 0;

St_S = 1/8;
Sv_S = 2/3;
lt_cm = 4;

rho = 1.225; % density [kg/m^3]

Nx=5; Ny=20;

% Assumptions:
% AR of wing and HTP is high enough (>6) -> lifting line can be assumed
% Ct/Cr (tip to chord) ratio is left free.

% Wing
lambda = 0.3;
A_ratio = 26;
cr_W=1; ct_W=lambda*cr_W; b_W=A_ratio*0.5*(cr_W+ct_W);
sweep_W=0; dihedral_W=0;
MGC=0.5*(cr_W+ct_W);

% %Horizontal tail
lambdah = lambda;
cr_H=0.5*cr_W; ct_H=lambda*cr_H; b_H=0.25*b_W;
sweep_H=0; dihedral_H=0; twist_H=0;

% %Vertical tail
cr_V=1*cr_H; ct_V=lambda*cr_V; b_V=4/3*b_H;
sweep_V=0; dihedral_V=0; twist_V=0; 

x_offset_W=0; z_offset_W=0;

%% Part 1: Compute ZL angle of wing for twist (0 to 8 deg) and CD.
% part1;

%% Part 2: Plotting wing's aerodynamic polar for alpha 0 to 10 deg.
% part2;

%% Part 3: Assumption -> Ground Effect. 
% Plot CL and CD for alpha 6deg, against AR 0.075Ao to 1.25Ao
% Ao is the nominal aspect ratio
%part3;

%% Part 4: CL and CD and Xcm for Sum(M)=0 when alpha = 6deg
part4;

%% Part 5: Assumption -> ground effect. CL, CD and CM_cm
part5;