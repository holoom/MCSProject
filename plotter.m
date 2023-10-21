%% plotter 

path = 'figs' ; 
if ~exist(ful, 'dir')
    mkdir(pwd, 'figs');
end
% Simulation Scenario 1

% scinario parameters
sim_num = 0; 
labels = {'under fixed SE current' , 'under varying current', 'with reference model', 'without reference model', 'following set of points'} ;
for i= 1:5
    close all
    sim_num = sim_num + 1; 

    if i == 1 
        simulation = 1 ; 
        ref_model = 1 ; 
        lbl = labels{i};
    elseif i== 2
        simulation = 2 ; 
        ref_model = 1 ;
        lbl = labels{i};
    elseif i ==3 
        simulation = 3 ; 
        ref_model = 1 ;
        lbl = labels{i};
    elseif i== 4
        simulation = 3 ; 
        ref_model = 0 ;
        lbl = labels{i};
    else 
        simulation = 4 ; 
        ref_model = 1 ;
        lbl = labels{i};
    end





% run the simulation 
sim('part1.slx') ; 

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

figure(1) 

plot(time, north)
hold on 
plot(time, east)
hold on 
plot(time, des_north)
hold on 
plot(time, des_east)
xlabel('Time (sec)')
ylabel('Amplitude (m)')
title(sprintf('DP %s (Translational motion)', lbl)); 
legend('North', 'East')
saveas(gcf, fullfile(path , sprintf('Translational_%02d.png', sim_num))); 

figure(2)

plot(time, psi) 
hold on 
plot(time, des_psi ) 
xlabel('Time (sec)')
ylabel('Amplitude (rad) ')
title(sprintf('DP %s (Rotational motion)', lbl)); 
saveas(gcf, fullfile(path , sprintf('Rotational_%02d.png', sim_num)));

figure(3) 


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

%scatter(selected_x, selected_y, 'filled');
%hold on 
quiver(selected_x, selected_y, selected_u, selected_v, 0.5, 'Color', 'r', 'LineWidth', 1.5); 


title('XY plot') 
xlabel('Amplitude in East direction (m)')
ylabel('Amplitude in North direction (m)')
legend('actual trajectory', 'desired trajectory')
saveas(gcf, fullfile(path , sprintf('xy_%02d.png', sim_num)));


 
end