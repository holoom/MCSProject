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
load('wind_coeff.mat')
% Initial position x, y, z, phi, theta, psi
eta0 = [0,0,0,0,0,0]';
% Initial velocity u, v, w, p, q, r
nu0 = [0,0,0,0,0,0]';
set_points_time = [1 2 3 4 5] * 1000 ; 

%% get omega and zeta from A matrix
A = vesselABC.Ar;
%% reference model parameters
t1 = 100;
t2 = 100; 
t3 = 100;
Af = diag ([1/t1, 1/t2, 1/t3]);

%
omega = [2*pi/600 2*pi/600 2*pi/600];
zeta = [0.9 0.9 0.9];
Gamma_Matrix = diag(omega.^2);
Omega_Matrix = diag (2*zeta.*omega);


%% Low-frequency control plant model
n_states = 15;
% mass matrix (rigid body + added mass)
M = [7.0101e6 0 0;
    0 8.5190e6 4.7187e5;
    0 4.7187e5 3.7973e9];
M_inv = inv(M);

% Damping matrix
D = [2.6486e5 0 0;
    0 8.8164e5 0;
    0 0 3.3774e8];
%% Wave-frequency control plant model
T_i = [9, 9, 9]; %s Ti, corresponding to wave frequency ωi = 2π/Ti, are in the range of 5 to 20 seconds in the North Sea for wind generated seas.
omega_i = 2*pi./T_i;


Sigma_matrix = diag([0.05, 0.05, 0.05]); % zeta ranging 0.05-0.1
Omega_matrix = diag(omega_i);

C_w = [zeros(3,3) eye(3)];
A_w = [zeros(3,3) eye(3);
       -Omega_matrix.^2 -2*Sigma_matrix*Omega_matrix];
Kw1 = 1;
Kw2 = 1;
Kw3 = 1;
K_w = diag ([Kw1 Kw2 Kw3]);
E_w = [zeros(3,3); K_w]; 

%% bias model 
T_b =60* diag([1 1 1]); % diagonal matrix based of the bais time constants
T_b_inv = inv(T_b);
E_b = diag ([1 1 1]); %is a diagonal scaling matrix
%% 
% Create a structure to store the matrices
model.n_states = n_states;
model.M_inv = M_inv;
model.M = M;
model.D = D;
model.A_w = A_w;
model.E_w = E_w;
model.C_w = C_w;
model.T_b = T_b;
model.E_b = E_b;

%% Environmental loads
% wind parameters 
mean_wind = 10 ; 
mean_wind_direction = 180 ; 
N =100;
z = 3;
n = 0.468;
U_10 = 12 ;
params_wind = [N z n U_10];

% current parameters
mean_current_speed = 0.2 ; 
mean_current_angle = 270 ; 

% waves parameters 
H_s = 2.5 ; 
T_p = 9   ; 
dir = 225 * pi/180 ;  % fromm north east 

%% Controller paramters
%surge, sway, yaw
Kp=[110700, 134500, 59964412];
Ki=[2213, 2690, 1199288];
Kd=[1233289, 1498741, 668054281];
% pr = 0.5 ;                           % 0.3 for pseudo
% Kp = [11400 78800 7830000*10]/pr; 
% Ki = [158.6 2300 106800*5]/pr;
% Kd = [203000 595000 76350000*5]/pr;




%% Observer paramters
%EKF design
x0 = zeros(15,1);
P0 = diag([1,1,1,pi/180,pi/180,pi/180,1,1,pi/180,1e5,1e5,1e5,1,1,pi/180]);

%Tunning Q and R
% Q =[100 0 0 0 0 0;          % best for 4 corner + quad (so far)/ not really good
%     0 100 0 0 0 0;
%     0 0 1 0 0 0;
%     0 0 0 0.1 0 0;
%     0 0 0 0 0.1 0;
%     0 0 0 0 0 0.00001] *0.001 ; % E[w' w]
% R =[0.3 0 0;
%     0 0.7 0;
%     0 0 0.0001]   ;


% Q = diag([0.4,0.4,0.01*pi/180,5e4,1e5,7e5].^2);
Q = diag([1,1,0.01*pi/180,5e4,5e4,7e5].^2);
%R = diag([0.5,0.1,3*pi/180].^2);
R = diag([5,5,3*pi/180].^2);



%Non linear observer
%choosing K1, K2, K3, and K4
T_peak = 25; % wave peak time period 
omega_wave_i = 2*pi/T_peak; % wave peak frequency
omega_wave_c = 1.3 * omega_wave_i;
zeta_ni = 1;
zeta_n = 0.1;
% K1 paramters 
k_1_3 = -2 * (zeta_ni - zeta_n) * (omega_wave_c/omega_wave_i);
k1 = k_1_3; k2 = k_1_3; k3 = k_1_3;
k_4_6 = 2 * omega_wave_i * (zeta_ni - zeta_n);
k4 = k_4_6; k5 = k_4_6; k6 = k_4_6;
% K2 paramters 
k_7_9 = omega_wave_c;
k7 = k_7_9; k8 = k_7_9; k9 = k_7_9;
% choose K4 in order of magnitude of Mass materix 
k13 = 7.0101e6; k14 = 8.5190e6; k15 = 3.7973e9;
% choose K3 less than 0.1*K4. K3 = 0.01 * K4 
%k13 = 1; k14 = 1; k15 = 1;
 
K1 =[diag([k1, k2, k3]);
    diag([k4, k5, k6])];
K2 = diag([k7, k8, k9]);
K4 = diag([k13, k14, k15]);
K3 = 0.01 * K4;
%% thruster allocation
xi = [thrusters.xposition];
yi = [thrusters.yposition];
thrust_rate = [thrusters.rate];
thrust_max = [thrusters.thrust];

%% Simulation Scenarios part 1
% sample_time = 0.01;
% simulation = 4; 
% ref_model = 1;
% observer = 1; % withoutobserver = 1, EKF = 2, NonlinearObserver = 3
%% Simulation Scenarios part 2
sample_time = 0.1;
set_points_time = [1 2 3 4 5] * 1000 ; 

% observer = 1; % withoutobserver = 1, EKF = 2, NonlinearObserver = 3
% 
% % enable params 
% observer = 1;
% ref_model = 1;    % 0 for 0 desired setpoint , 1 for 4 corners     
% use_ref = 1;     % 0 for desired setpoint without reference model , 1 use reference model ()
% enable_waves = 0;
% use_thrust_allocation  = 1 ;            % 0 without thrust allocation , 1 with thrust allocation 
% use_env_forces = 0 ;      % 0 for no environmental forces , 1 for enalbling the environmental forces
% allocation_method = 1 ;    % 0 for quadratic programming method , and 1 for pseudo inverse method
% desired_DP_force = [1 1 1]*1e4;
% use_fixed_con = 0 ; % 1 for fixed desired DP force, 0 for using PID controller



% Simulation 4 - Observer selection
% enable params
observer = 2; % withoutobserver = 1, EKF = 2, NonlinearObserver = 3
ref_model = 1;    % 0 for 0 desired setpoint , 1 for 4 corners     
use_ref = 1 ;     % 0 for desired setpoint without reference model , 1 use reference model ()
use_thrust_allocation  = 1 ;            % 0 without thrust allocation , 1 with thrust allocation 
use_env_forces = 1;      % 0 for no environmental forces , 1 for enalbling the environmental forces
enable_waves = 1;        % 0 for not using waves, 1 for enabling the waves
allocation_method = 0;   % 0 for quadratic programming method , and 1 for pseudo inverse method
desired_DP_force = [1 1 1]*1e4;
use_fixed_con = 0 ; % 1 for fixed desired DP force, 0 for using PID controller


%% Simulation Scenarios and results 

path = 'figs' ; 
if ~exist(fullfile(pwd, path), 'dir')
    mkdir(pwd, 'figs');
end

