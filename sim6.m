path = 'figs' ; 
if ~exist(path, 'dir')
    mkdir(pwd, path);
end 


sub_path = 'figs/sim6' ; 
if ~exist(sub_path, 'dir')
    mkdir(pwd, sub_path);
end 



% sim4 parameters
sim_time= 200 ;


observer = 1; % withoutobserver = 1, EKF = 2, NonlinearObserver = 3
use_env_forces = 1;      % 0 for no environmental forces , 1 for enalbling the environmental forces
enable_waves = 1;        % 0 for not using waves, 1 for enabling the waves
use_fixed_con = 0 ; % 1 for fixed desired DP force, 0 for using PID controller
desired_DP_force = [1 1 1]*10^4;
% unnecessary parameters
use_thrust_allocation  = 1 ;            % 0 without thrust allocation , 1 with thrust allocation 
allocation_method = 0;   % 0 for quadratic programming method , and 1 for pseudo inverse method
ref_model = 0;    % 0 for 0 desired setpoint , 1 for 4 corners     
use_ref = 1 ;     % 0 for desired setpoint without reference model , 1 use reference model ()

% environmental params 
% wind parameters 
mean_wind = 10 ; 
N =100;
z = 3;
n = 0.468;
U_10 = 12 ;
params_wind = [N z n U_10];

% current parameters
mean_current_speed = 0.2 ; 

% waves parameters 
H_s = 4 ; 
T_p = 8   ; 

% Define the maximum thrust values for each thruster
thrust_max = [125000, 150000, 125000, 300000, 300000];

% Pre-define the environmental directions for which to run the simulations
directions = 0:10:350; % Directions from 0 to 350 degrees with 10-degree increments

% Initialize arrays to store the average thrust utilization
avg_thrust_utilization = zeros(size(directions));
avg_thrust_utilization_with_constraints = zeros(size(directions));

% Loop through each direction and run the simulation
for i = 1:length(directions)
    % Change direction
    mean_current_angle = directions(i); % Update the mean current angle
    mean_wind_direction = directions(i); % Update the mean wind direction
    dir = directions(i) * pi/180; % Convert direction to radians

    sim('part1.slx');

    % Get the thrust output and eta from the simulation
    out_thrust = logsout.getElement('out_thrust') ; 
    out_thrust = out_thrust.Values.Data;
    eta = logsout.getElement('eta') ;
    north = eta.Values.Data(: , 1) ; 
    east  = eta.Values.Data(: , 2) ;
    psi   = eta.Values.Data(: , 3) ; 

    % Calculate the average thrust utilization for each thruster
    for j = 1:size(out_thrust, 2) % Assuming thrust_out is a matrix with columns corresponding to each thruster
        avg_thrust_utilization(i) = avg_thrust_utilization(i) + mean(out_thrust(:, j)) / thrust_max(j);
    end
    avg_thrust_utilization(i) = (avg_thrust_utilization(i)/5)*100;
    % Check if the position and heading deviations are within the constraints
    position_deviation = max(sqrt(north.^2 + east.^2));
    heading_deviation = max(abs(psi)); % Assuming psi is in radians

    if position_deviation <= 3 && heading_deviation <= (3 * pi/180)
        avg_thrust_utilization_with_constraints(i) = avg_thrust_utilization(i);
    else
        avg_thrust_utilization_with_constraints(i) = NaN; % Not displaying values that exceed constraints
    end
end

% Plot the average thrust utilization without constraints
figure;
polarplot(deg2rad(directions), avg_thrust_utilization,'bo-');
title('Average Thrust Utilization Without Constraints [%]');
% Saving the figure
saveas(gcf, fullfile(sub_path, 'Average Thrust Utilization Without Constraint.png'));
% Plot the average thrust utilization with constraints
figure;
polarplot(deg2rad(directions), avg_thrust_utilization_with_constraints, 'ro-');
title('Average Thrust Utilization With Constraints [%]');
% Saving the figure
saveas(gcf, fullfile(sub_path, 'Average Thrust Utilization With Constraint.png'));
