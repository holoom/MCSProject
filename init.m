%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init()                                                                  %
%                                                                         %              
% Set initial parameters for part1.slx and part2.slx                      %
%                                                                         %
% Created:      2018.07.12	Jon Bjørnø                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

load('supply.mat');
load('supplyABC.mat');
load('thrusters_sup.mat')

% Initial position x, y, z, phi, theta, psi
eta0 = [0,0,0,0,0,0]';
% Initial velocity u, v, w, p, q, r
nu0 = [0,0,0,0,0,0]';

%% get omega and zeta from A matrix
A = vesselABC.Ar;
%% reference model parameters
t1 = 10;
t2 = 10; 
t3 = 20;
Af = diag ([1/t1, 1/t2, 1/t3]);

%
omega = [2*pi/100 2*pi/100 2*pi/50];
zeta = [0.9 0.7 1];
Gamma_Matrix = diag(omega.^2);
Omega_Matrix = diag (2*zeta.*omega);
%% Controller paramters
%surge, sway, yaw
Kp = [1 1 1]; 
Kd = [1 1 1];
Ki = [1 1 1];

%% PID tunning usign LQG


