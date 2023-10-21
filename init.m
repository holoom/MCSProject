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
omega = [2*pi/100 2*pi/100 2*pi/100];
zeta = [0.9 0.9 0.9];
Gamma_Matrix = diag(omega.^2);
Omega_Matrix = diag (2*zeta.*omega);
%% Controller paramters
%surge, sway, yaw
Kp = [114000 788000 78300000]; 
Ki = [1586 23000 1068000];
Kd = [2030000 5950000 763500000];


%% PID tunning usign LQG


%% Simulation Scenarios
sample_time = 0.01;
simulation = 4; 
ref_model = 1;

%% Simulation Scenarios and results 

path = 'figs' ; 
if ~exist(fullfile(pwd, path), 'dir')
    mkdir(pwd, 'figs');
end

