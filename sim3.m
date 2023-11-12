%%%%%%%%% this is simulation 3 %%%%%%%%%


path = 'figs' ; 
if ~exist(path, 'dir')
    mkdir(pwd, path);
end 


sub_path = 'figs/sim3' ; 
if ~exist(sub_path, 'dir')
    mkdir(pwd, sub_path);
end 



observer = 1;                 % no observer
use_fixed_con = 0 ;           % use DP controller  
use_env_forces = 1 ; 
use_waves = 1 ; 
allocation_method = 1 ;       % pseudo inverse
ref_model = 1 ; 
use_ref = 1 ; 


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
title('DP with unfiltered environmental loads (Translational motion)'); 
legend('North', 'East')
saveas(gcf, fullfile(sub_path  , 'Translational_sim3.png')); 

figure

plot(time, psi) 
hold on 
plot(time, des_psi ) 
xlabel('Time (sec)')
ylabel('Amplitude (rad) ')
title('DP with unfiltered environmental loads (Rotational motion)'); 
saveas(gcf, fullfile(sub_path  , 'Rotational_sim3.png')); 

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
saveas(gcf, fullfile(sub_path  , 'xy_sim3.png')); 





