% Define the parameters and variables (replace with actual values)
m = 1.0;      % Placeholder for m

Iz = 1.0;     % Placeholder for Iz

Xu_dot = 0;   % Placeholder for Xu˙
Yv_dot = 0;   % Placeholder for Yv˙
Yr_dot = 0;   % Placeholder for Yr˙
Nv_dot = 0;   % Placeholder for Nv˙
Nr_dot = 0;   % Placeholder for Nr˙

xG = 0;      % Placeholder for xG

Xu = 0;       % Placeholder for -Xu
Yv = 0;       % Placeholder for -Yv
Yr = 0;       % Placeholder for -Yr
Nv = 0;       % Placeholder for -Nv
Nr = 0;       % Placeholder for -Nr

%% 
% Define the matrices
M = [m-Xu_dot 0 0;
     0 m-Yv_dot m*xG-Yr_dot;
     0 m*xG-Nv_dot Iz-Nr_dot];
M_inv = inv(M);
D = [-Xu 0 0;
     0 -Yv -Yr;
     0 -Nv -Nr];
%% 
sigma = diag([1, 1, 1]);
omega = diag([1, 1, 1]);

C_w = [zeros(3,3) eye(3)];
A_w = [zeros(3,3) eye(3);
       -omega.^2 -2*sigma*omega];
Kw1 = 1;
Kw2 = 1;
Kw3 = 1;
K_w = diag ([Kw1 Kw2 Kw3]);
E_w = [zeros(3,3); K_w]; 

%% bias
T_b = diag([1 1 1]); % diagonal matrix based of the bais time constants
E_b = diag ([1 1 1]); %is a diagonal scaling matrix


%% EKF design
n_states = 15;
sample_time = 0.1;


x0 = zeros(15,1);
P0 = eye(15);

%Tunning Q and R
Q = eye(6); % E[w' w]
R = eye(3); % E[v' v]
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









