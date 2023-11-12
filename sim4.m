path = 'figs' ; 
if ~exist(path, 'dir')
    mkdir(pwd, path);
end 


sub_path = 'figs/sim4' ; 
if ~exist(sub_path, 'dir')
    mkdir(pwd, sub_path);
end 



% sim4 parameters
sim_time= 200 ;
observer = 1; % withoutobserver = 1, EKF = 2, NonlinearObserver = 3
use_env_forces = 1;      % 0 for no environmental forces , 1 for enalbling the environmental forces
enable_waves = 1;        % 0 for not using waves, 1 for enabling the waves
use_fixed_con = 1 ; % 1 for fixed desired DP force, 0 for using PID controller
desired_DP_force = [1 1 1]*10^4;
% unnecessary parameters
use_thrust_allocation  = 0 ;            % 0 without thrust allocation , 1 with thrust allocation 
allocation_method = 1;   % 0 for quadratic programming method , and 1 for pseudo inverse method
ref_model = 1;    % 0 for 0 desired setpoint , 1 for 4 corners     
use_ref = 1 ;     % 0 for desired setpoint without reference model , 1 use reference model ()

% environmental params 
% wind parameters 
mean_wind = 10 ; 
mean_wind_direction = 180 ;  %  from north
N =100;
z = 3;
n = 0.468;
U_10 = 12 ;
params_wind = [N z n U_10];

% current parameters
mean_current_speed = 0.2 ; 
mean_current_angle = 270 ;     % from east

% waves parameters 
H_s = 2.5 ; 
T_p = 9   ; 
dir = 225 * pi/180 ;  % fromm north east 

enable_waves = 1;        % 0 for not using waves, 1 for enabling the waves
% run the simulation 
sim('part1.slx') ; 


% plotting


% data extraction 
eta = logsout.getElement('eta') ;
eta_ekf = logsout.getElement('eta_ekf') ;
eta_nlo = logsout.getElement('eta_nlo') ;

nu = logsout.getElement('nu') ;
nu_ekf = logsout.getElement('nu_ekf') ;
nu_nlo = logsout.getElement('nu_nlo') ;
time= eta.Values.Time ; 

% plotting

figure;

% Subplot 1: North component
subplot(3, 1, 1);
plot(time, eta.Values.Data(:, 1));
hold on;
plot(time, eta_ekf.Values.Data(1, :));
plot(time, eta_nlo.Values.Data(:,1));
xlabel('Time (sec)');
ylabel('Amplitude (m)');
legend('Real', 'EKF', 'NPO');
title('North Component');

% Subplot 2: East component
subplot(3, 1, 2);
plot(time, eta.Values.Data(:, 2));
hold on;
plot(time, eta_ekf.Values.Data(2,:));
plot(time, eta_nlo.Values.Data(:,2));
xlabel('Time (sec)');
ylabel('Amplitude (m)');
legend('Real', 'EKF', 'NPO');
title('East Component');

% Subplot 3: Psi component
subplot(3, 1, 3);
plot(time, eta.Values.Data(:, 3));
hold on;
plot(time, eta_ekf.Values.Data(3,:));
hold on;
plot(time,eta_nlo.Values.Data(:,3));
xlabel('Time (sec)');
ylabel('Amplitude (rad)');
legend('Real', 'EKF', 'NPO');
title('Psi Component');

% Setting a title for the entire figure
sgtitle('Observer eta Estimation (Including Waves)');

% Saving the figure
saveas(gcf, fullfile(sub_path, 'eta_sim4.png')); 

figure;

% Subplot 1: North component
subplot(3, 1, 1);
plot(time, nu.Values.Data(1,:));
hold on;
plot(time, nu_ekf.Values.Data(1, :));
plot(time, nu_nlo.Values.Data(:, 1));
xlabel('Time (sec)');
ylabel('Amplitude (m/s)');
legend('Real', 'EKF', 'NPO');
title('North Component');

% Subplot 2: East component
subplot(3, 1, 2);
plot(time, nu.Values.Data(2,:));
hold on;
plot(time, nu_ekf.Values.Data(2,:));
plot(time, nu_nlo.Values.Data(:, 2));
xlabel('Time (sec)');
ylabel('Amplitude (m/s)');
legend('Real', 'EKF', 'NPO');
title('East Component');

% Subplot 3: Psi component
subplot(3, 1, 3);
plot(time, nu.Values.Data(3,:));
hold on;
plot(time, nu_ekf.Values.Data(3,:));
plot(time, nu_nlo.Values.Data(:, 3));
xlabel('Time (sec)');
ylabel('Amplitude (rad/s)');
legend('Real', 'EKF', 'NPO');
title('Psi Component');

% Setting a title for the entire figure
sgtitle('Velocity Estimation (Including Waves)');

% Saving the figure
saveas(gcf, fullfile(sub_path, 'nu_sim4.png')); 

% xyplot
figure;
plot(eta.Values.Data(:, 1), eta.Values.Data(:, 2));
hold on;
plot(eta_ekf.Values.Data(1,:), eta_ekf.Values.Data(2,:));
hold on;
plot(eta_nlo.Values.Data(:, 1), eta_nlo.Values.Data(:, 2));
xlabel('North (m)');
ylabel('East (m)');
legend('Real', 'EKF', 'NPO');
title('XY-plot for Observer Output (Including Waves)');
% Saving the figure
saveas(gcf, fullfile(sub_path, 'xy_sim4.png'));

enable_waves = 0;        % 0 for not using waves, 1 for enabling the waves
% run the simulation 
sim('part1.slx') ; 


% plotting


% data extraction 
eta = logsout.getElement('eta') ;
eta_ekf = logsout.getElement('eta_ekf') ;
eta_nlo = logsout.getElement('eta_nlo') ;

nu = logsout.getElement('nu') ;
nu_ekf = logsout.getElement('nu_ekf') ;
nu_nlo = logsout.getElement('nu_nlo') ;
time= eta.Values.Time ; 

% plotting

figure;

% Subplot 1: North component
subplot(3, 1, 1);
plot(time, eta.Values.Data(:, 1));
hold on;
plot(time, eta_ekf.Values.Data(1, :));
plot(time, eta_nlo.Values.Data(:,1));
xlabel('Time (sec)');
ylabel('Amplitude (m)');
legend('Real', 'EKF', 'NPO');
title('North Component');

% Subplot 2: East component
subplot(3, 1, 2);
plot(time, eta.Values.Data(:, 2));
hold on;
plot(time, eta_ekf.Values.Data(2,:));
plot(time, eta_nlo.Values.Data(:,2));
xlabel('Time (sec)');
ylabel('Amplitude (m)');
legend('Real', 'EKF', 'NPO');
title('East Component');

% Subplot 3: Psi component
subplot(3, 1, 3);
plot(time, eta.Values.Data(:, 3));
hold on;
plot(time, eta_ekf.Values.Data(3,:));
hold on;
plot(time,eta_nlo.Values.Data(:,3));
xlabel('Time (sec)');
ylabel('Amplitude (rad)');
legend('Real', 'EKF', 'NPO');
title('Psi Component');

% Setting a title for the entire figure
sgtitle('Observer eta Estimation (Without Waves)');

% Saving the figure
saveas(gcf, fullfile(sub_path, 'eta_sim4_without_waves.png')); 

figure;

% Subplot 1: North component
subplot(3, 1, 1);
plot(time, nu.Values.Data(1,:));
hold on;
plot(time, nu_ekf.Values.Data(1, :));
plot(time, nu_nlo.Values.Data(:, 1));
xlabel('Time (sec)');
ylabel('Amplitude (m/s)');
legend('Real', 'EKF', 'NPO');
title('North Component');

% Subplot 2: East component
subplot(3, 1, 2);
plot(time, nu.Values.Data(2,:));
hold on;
plot(time, nu_ekf.Values.Data(2,:));
plot(time, nu_nlo.Values.Data(:, 2));
xlabel('Time (sec)');
ylabel('Amplitude (m/s)');
legend('Real', 'EKF', 'NPO');
title('East Component');

% Subplot 3: Psi component
subplot(3, 1, 3);
plot(time, nu.Values.Data(3,:));
hold on;
plot(time, nu_ekf.Values.Data(3,:));
plot(time, nu_nlo.Values.Data(:, 3));
xlabel('Time (sec)');
ylabel('Amplitude (rad/s)');
legend('Real', 'EKF', 'NPO');
title('Psi Component');

% Setting a title for the entire figure
sgtitle('Velocity Estimation (Without Waves)');

% Saving the figure
saveas(gcf, fullfile(sub_path, 'nu_sim4_without_waves.png')); 

% xyplot
figure;
plot(eta.Values.Data(:, 1), eta.Values.Data(:, 2));
hold on;
plot(eta_ekf.Values.Data(1,:), eta_ekf.Values.Data(2,:));
hold on;
plot(eta_nlo.Values.Data(:, 1), eta_nlo.Values.Data(:, 2));
xlabel('North (m)');
ylabel('East (m)');
legend('Real', 'EKF', 'NPO');
title('XY-plot for Observer Output (Without Waves)');
% Saving the figure
saveas(gcf, fullfile(sub_path, 'xy_sim4_without_waves.png'));
