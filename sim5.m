%%%%%%%%% this is simulation 5 %%%%%%%%%


path = 'figs' ; 
if ~exist(path, 'dir')
    mkdir(pwd, path);
end 


sub_path = 'figs/sim5' ; 
if ~exist(sub_path, 'dir')
    mkdir(pwd, sub_path);
end 



% sim2 parameters
observer = 2; % withoutobserver = 1, EKF = 2, NonlinearObserver = 3
use_env_forces = 1;      % 0 for no environmental forces , 1 for enalbling the environmental forces
enable_waves = 1;        % 0 for not using waves, 1 for enabling the waves
allocation_method = 1;   % 0 for quadratic programming method , and 1 for pseudo inverse method
ref_model = 1;    % 0 for 0 desired setpoint , 1 for 4 corners     
use_ref = 1 ;     % 0 for desired setpoint without reference model , 1 use reference model ()
use_thrust_allocation  = 1 ;            % 0 without thrust allocation , 1 with thrust allocation 
use_fixed_con = 0 ; % 1 for fixed desired DP force, 0 for using PID controller

% unnecessary parameters
desired_DP_force = [0 0 0];


sim_time= 6000 ; 


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


% run the simulation 
sim('part1.slx') ; 


% plotting


% data extraction 
eta =          logsout.getElement('eta') ; 
desired_eta =  logsout.getElement('set_point') ;
time= eta.Values.Time ; 

% current pose
north = eta.Values.Data(: , 1) ; 
east  = eta.Values.Data(: , 2) ;
psi   = eta.Values.Data(: , 3) ; 

% desired pose 
des_north = desired_eta.Values.Data(: , 1) ;
des_east  = desired_eta.Values.Data(: , 2) ;
des_psi   = desired_eta.Values.Data(: , 3) ; 


% plotting

figure

plot(time, north)
hold on 
plot(time, east)
hold on 
plot(time, des_north)
hold on 
plot(time, des_east)
xlabel('Time (sec)')
ylabel('Amplitude (m)')
title('Full DP system (Translational motion)'); 
legend('North', 'East')
saveas(gcf, fullfile(sub_path  , 'Translational_sim5.png')); 

figure

plot(time, psi) 
hold on 
plot(time, des_psi ) 
xlabel('Time (sec)')
ylabel('Amplitude (rad) ')
title('Full DP system (Rotational motion)'); 
saveas(gcf, fullfile(sub_path  , 'Rotational_sim5.png')); 

figure


plot(east, north) ; 
hold on
plot(des_east, des_north), 
hold on
indices = 1:1000:length(east) ;

selected_x = east(indices);
selected_y = north(indices);
selected_u = cos(-psi + (pi/2)) ; 
selected_u = selected_u(indices);
selected_v = sin(-psi + (pi/2)) ; 
selected_v = selected_v(indices); 


quiver(selected_x, selected_y, selected_u, selected_v, 0.5, 'Color', 'r', 'LineWidth', 1.5); 


title('XY plot') 
xlabel('Amplitude in East direction (m)')
ylabel('Amplitude in North direction (m)')
legend('actual trajectory', 'desired trajectory')
saveas(gcf, fullfile(sub_path  , 'xy_sim5.png')); 







