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

%% reference model parameters
t1 = 2;
t2 = 8; 
t3 = 10;
Af = diag ([1/t1, 1/t2, 1/t3]);

%
omega = [2 3 4];
zeta = [0.1 0.2 0.3];
Gamma_Matrix = diag(omega.^2);
Omega_Matrix = diag (2*zeta.*omega);
%% Controller paramters
%surge, sway, yaw
Kp = [1 1 1]; 
Kd = [1 1 1];
Ki = [1 1 1];

%% PID tunning usign LQG


