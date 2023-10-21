syms x y psi u v r; 
syms b1 b2 b3;
syms zeta1 zeta2 zeta3 zeta4 zeta5 zeta6
syms A_w 
b = [b1 b2 b3]';
eta = [x y psi]';
nu = [u v r]';
zeta = [zeta1 zeta2 zeta3 zeta4 zeta5 zeta6]';

%% 
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
c = cos(psi);
s = sin(psi);

% Create the rotation matrix R(ψ)
R_psi = [c, -s, 0;
     s, c, 0;
     0, 0, 1];


% Define the matrices
M = [m-Xu_dot 0 0;
     0 m-Yv_dot m*xG-Yr_dot;
     0 m*xG-Nv_dot Iz-Nr_dot];

D = [-Xu 0 0;
     0 -Yv -Yr;
     0 -Nv -Nr];
%% EoM
tau = zeros(3,1);
eta_dot = R_psi*nu;
nu_dot = M\(tau + R_psi' * b - D * nu); 
%% 
sigma = diag([1, 1, 1]);
omega = diag([1, 1, 1]);
W_w = randn(6,1); % zero mean gaussian white noise  

C_w = [zeros(3,3) eye(3)];
A_w = [zeros(3,3) eye(3);
       -omega.^2 -2*sigma*omega];
Kw1 = 1;
Kw2 = 1;
Kw3 = 1;
K_w = diag ([Kw1 Kw2 Kw3]);
E_w = [zeros(3,3); K_w]; 
zeta_dot = A_w * zeta+ E_w * W_w;
eta_w = C_w * zeta;

%% bias
w_b = randn(3,1); % zeros mean gaussian noise 
T_b = diag([1 1 1]); % diagonal matrix based of the bais time constants
E_b = diag ([1 1 1]); %is a diagonal scaling matrix
b_dot = -inv(T_b)*b + E_b * w_b; %markov process model T_b is 3x3 and E_b is 3x3
%or
% b_dot = E_b * w_b; %wiener process model

%% Measurement model
v_m = zeros(3,1);
y = eta + C_w * zeta + v_m; %v_m is the measurement noise
%% 
f_x = [A_w * zeta;
    R * nu;
    - T_b/b;
    M\(R_psi' * b - D * nu)];



%% EKF design
n_states = 15;
%simulation time
Tsim = 1000;
dt = 0.1;
n_samples = Tsim/dt;


x = zeros(n_states, n_samples+1);
p = zeros(n_states, n_states, n_samples+1);
X = zeros(n_states, n_samples);
%intiallization 
x(:,1) = x0;
P(:,:,1) = P0; % E[(x(0) − x^(0))(x(0) − x^(0))']

%Tunning Q and R
Q = eye(6); % E[w' w]
R = eye(3); % E[v' v]

for i = 1:n_samples
P_k =P(:,:,i) ;
x_k = x(:,i);
y_k = zeros(3,1);

%Corrector
K_k = P_k * H' * (H * P_k * H')^(-1);
P_k_ = (eye(n_states) - K_k * H) * P_k * (eye(n_states) - K_k * H)' + K_k * R * K_k';
x_k_ = x_k + K_k * (y_k - H * x_k);
X(:,i) = x_k_;
%predictor
P(:,:,i+1) = Phi_k * P_k_ * Phi_k' + Gamma_k * Q * Gamma_k';
x(:,i+1) = f_k(x_k_, u_k);
end


%% 












