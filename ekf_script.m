
y_k = [1 1 1]';
tau_est = [1 1 1]';
%model
%x = [zeta eta bias nu]

%intiallization at 
x_k = x0;

P_k = P0;

M_inv = model.M_inv;
B = [zeros(6,3);
    zeros(3,3);
    zeros(3,3);
    M_inv];

E = [model.E_w zeros(6,3);
    zeros(3,3) zeros(3,3);
    zeros(3,3) model.E_b;
    zeros(3,3) zeros(3,3)];

H = [model.C_w eye(3) zeros(3,3) zeros(3,3)];


%Corrector
K_k = P_k * H' * (H * P_k * H')^(-1);
P_k_ = (eye(model.n_states) - K_k * H) * P_k * (eye(model.n_states) - K_k * H)' + K_k * R * K_k';
x_k_ = x_k + K_k * (y_k - H * x_k);


%%
zeta = x_k_(1:6); eta = x_k_(7:9); b = x_k_(10:12); nu = x_k_(13:15);
%rotation matrix
psi = eta(3);
c = cos(psi);
s = sin(psi);

% Create the rotation matrix R(ψ)
R_psi = [c, -s, 0;
     s, c, 0;
     0, 0, 1];

T_b_inv = inv(model.T_b);
f_x = [model.A_w * zeta;
        R_psi * nu;
        -T_b_inv * b;
        -M_inv * model.D * nu + M_inv * R_psi' * b];

jac_f_x = [model.A_w zeros(6,3) zeros(6,3) zeros(6,3);
           zeros(3,6) zeros(3,3) zeros(3,3) R_psi;
           zeros(3,6) zeros(3,3) -T_b_inv zeros(3,3);
           zeros(3,6) zeros(3,3) M_inv * R_psi' -M_inv * model.D];

%parial derivative of f_x with respect to psi
% jac_f_x(7,9) = - sin(psi) * u - cos(psi) * v ;
% jac_f_x(8,9) = cos(psi) * u - sin(psi) * v ;
% jac_f_x(9,9) = 0 ;
diff_R = [-s -c 0;
          c -s 0;
          0 0 0];
jac_f_x(7:9,9) = diff_R * nu;
jac_f_x(13:15,9) = M_inv * diff_R' * b ;


Phi_k = eye(model.n_states) + sample_time * jac_f_x;
Gamma_k = sample_time * E;

%%
%predictor
P_k = Phi_k * P_k_ * Phi_k' + Gamma_k * Q * Gamma_k';
x_k = x_k_ + sample_time * (f_x + B * tau_est );


eta_hat = x_k_(7:9) ;
nu_hat = x_k_(13:15);
