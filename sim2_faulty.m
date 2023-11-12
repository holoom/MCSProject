%%%%%%%%%% sim 2 with thruster 2 and 4 are disabled and usin pseudo
%%%%%%%%%% inverse allocator
 



sub_path = 'figs/sim2/faulty_thrusters' ; 
if ~exist(sub_path, 'dir')
    mkdir(pwd, sub_path);
end 
% sim2 parameters
observer = 1; % withoutobserver = 1, EKF = 2, NonlinearObserver = 3
use_env_forces = 0;      % 0 for no environmental forces , 1 for enalbling the environmental forces
enable_waves = 0;        % 0 for not using waves, 1 for enabling the waves
allocation_method = 1;   % 0 for quadratic programming method , and 1 for pseudo inverse method
ref_model = 1;    % 0 for 0 desired setpoint , 1 for 4 corners     
use_ref = 1 ;     % 0 for desired setpoint without reference model , 1 use reference model ()
use_thrust_allocation  = 1 ;            % 0 without thrust allocation , 1 with thrust allocation 
use_fixed_con = 0 ; % 1 for fixed desired DP force, 0 for using PID controller

% unnecessary parameters
desired_DP_force = [0 0 0];
sim_time= 6000 ; 


% simulate 
sim('part1.slx') 

% extract information

eta = logsout.getElement('eta') ; 
desired_eta =  logsout.getElement('set_point');
desired_eta = desired_eta.Values.Data ;
tau_d = logsout.getElement('tau_d');
tau_d = tau_d.Values.Data ; 
out_thrust = logsout.getElement('out_thrust') ; 
out_thrust = out_thrust.Values.Data; 
alphas = logsout.getElement('alphas'); 
alphas = alphas.Values.Data; 
thrusters_dynamics = logsout.getElement('thrusters_dynamics');
thrusters_dynamics = thrusters_dynamics.Values.Data; 
des_nu = logsout.getElement('nu_d');
des_nu = des_nu.Values.Data; 
nu = logsout.getElement('nu');
nu = nu.Values.Data;
nu = reshape(nu, [size(des_nu)]) ; 
time= eta.Values.Time ; 

% current pose
north = eta.Values.Data(: , 1) ; 
east  = eta.Values.Data(: , 2) ;
psi   = eta.Values.Data(: , 3) * 180/pi ;  


l = {'Surge' , 'Sway' , 'heave' ,'roll' ,'pitch' , 'Yaw'  } ; 

for i = 1:min(size(tau_d))
    f = tau_d(: , i) ; 
    ff = thrusters_dynamics(: , i) ;
    
    figure 

    plot(time , f)
    hold on
    plot(time, ff) 

    xlabel('Time (sec)')
    ylabel('Magnitude (N) ')
    legend('desired thrust' , 'allocated thrust')
    title(sprintf('Performance of thrust allocation in %s', l{i})); 
    saveas(gcf, fullfile(sub_path , strcat( num2str(i), sprintf('fpallocation_sim%02d.png', 2))));
end 


for i = 1:min(size(alphas))
    f = alphas(: , i) ; 
    ff = out_thrust(: , i) ;
    
    figure 

    plot(time , f)
    xlabel('Time (sec)')
    ylabel('Angle (rad) ')
    title(sprintf(' Direction of thrust in thruster %0d', i)); 
    
    saveas(gcf, fullfile(sub_path , strcat( num2str(i) , sprintf('fpoutput_thrust_direction_sim_%02d.png', 2))));
    
    figure
    
    plot(time, ff) 

    xlabel('Time (sec)')
    ylabel('Thrust (N) ')
    title(sprintf('Output thrust by thruster %02d', i)); 
    saveas(gcf, fullfile(sub_path , strcat( num2str(i) , sprintf('fpoutput_thrust_magnitude_sim_%02d.png', 2))));
end 

plot(desired_eta(: , 2) , desired_eta(: , 1))
hold on
plot(east, north) ; 

hold on
indices = 1:1000:length(east) ;

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
legend('desired trajectory' , 'ship trajectory')
saveas(gcf, fullfile(sub_path , sprintf('pxy_faulty_thrusters_%02d.png', 2))); 