
path = 'figs' ; 
if ~exist(path, 'dir')
    mkdir(pwd, path);
end 


sub_path = 'figs/sim1' ; 
if ~exist(sub_path, 'dir')
    mkdir(pwd, sub_path);
end 


% sim1 parameters

observer = 1; % withoutobserver = 1, EKF = 2, NonlinearObserver = 3
use_env_forces = 1;      % 0 for no environmental forces , 1 for enalbling the environmental forces
enable_waves = 1;        % 0 for not using waves, 1 for enabling the waves
allocation_method = 1;   % 0 for quadratic programming method , and 1 for pseudo inverse method
desired_DP_force = [0 0 0];
use_fixed_con = 1 ; % 1 for fixed desired DP force, 0 for using PID controller
thrusterfault = 1 ;    % healthy 

% unnecessary parameters
ref_model = 1;    % 0 for 0 desired setpoint , 1 for 4 corners     
use_ref = 1 ;     % 0 for desired setpoint without reference model , 1 use reference model ()
use_thrust_allocation  = 1 ;            % 0 without thrust allocation , 1 with thrust allocation 



sim_time= 300 ;

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

% data extraction 
eta = logsout.getElement('eta') ; 
time= eta.Values.Time ; 

% current pose
north = eta.Values.Data(: , 1) ; 
east  = eta.Values.Data(: , 2) ;
psi   = eta.Values.Data(: , 3) ;  

% plotting

figure(1) 

plot(time, north)
hold on 
plot(time, east)
hold on
plot(time, psi) 
xlabel('Time (sec)')
ylabel('Distance (m)')
title('Drifing Vessel under Environmental Loads (Position)'); 
legend('North', 'East')
saveas(gcf, fullfile(sub_path , sprintf('Translational_%02d.png', 1))); 

figure(2)

plot(time, psi ) 

xlabel('Time (sec)')
ylabel('Amplitude (rad) ')
title('Drifing Vessel under Environmental Loads (Heading)'); 
saveas(gcf, fullfile(sub_path , sprintf('Rotational_%02d.png', 1)));


figure(3) 

plot(east, north) ; 
hold on
indices = 1:500:length(east) ;

selected_x = east(indices);
selected_y = north(indices);
selected_u = cos(-psi + (pi/2)) ; 
selected_u = selected_u(indices);
selected_v = sin(-psi + (pi/2)) ; 
selected_v = selected_v(indices); 

quiver(selected_x, selected_y, selected_u, selected_v, 0.2, 'Color', 'r', 'LineWidth', 1.5); 


title('XY plot') 
xlabel('Position in East direction (m)')
ylabel('Position in North direction (m)')
legend('ship path')
saveas(gcf, fullfile(sub_path , sprintf('xy_%02d.png', 1)));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



% extra tests 
use_waves= 0 ; 
mean_current_speed = 0  ; 
mean_wind = 40 ; 
mean_wind_direction = 135 ; 
sim_time = 2000;  

sim('part1.slx') 
 
% state 

eta = logsout.getElement('eta') ; 
time= eta.Values.Time ; 

% current pose
north = eta.Values.Data(: , 1) ; 
east  = eta.Values.Data(: , 2) ;
psi   = eta.Values.Data(: , 3) * 180/pi ;  

figure(4)

plot(time, psi ) 

xlabel('Time (sec)')
ylabel('Amplitude (degree) ')
title('Heading of the Vessel under strong wind in South East'); 
saveas(gcf, fullfile(sub_path , strcat( num2str(1) , sprintf('extra%02d.png', 1)))); 

use_waves= 0 ; 
mean_current_speed = 40  ;
mean_current_angle = 135  ; 
mean_wind = 0 ; 
sim_time = 2000; 


sim('part1.slx') 
 
% state 

eta = logsout.getElement('eta') ; 
time= eta.Values.Time ; 

% current pose
north = eta.Values.Data(: , 1) ; 
east  = eta.Values.Data(: , 2) ;
psi   = eta.Values.Data(: , 3) * 180/pi ;  

figure(5)

plot(time, psi ) 

xlabel('Time (sec)')
ylabel('Amplitude (degree) ')
title('Heading of the Vessel under strong current in South East'); 
saveas(gcf, fullfile(sub_path , strcat( num2str(2) , sprintf('extra%02d.png', 1))));

