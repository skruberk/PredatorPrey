%muton and glutton 
clc;    % Clear the command window.
clearvars;
close all;  % Close all figs
workspace;  % Make sure the workspace panel is showing.
% Parameters
num_glutt = 5; % Number of Gluttons
grid_size = [100, 100]; % Size of the 2D world
initial_girth = 5; % Initial size of Gluttons
metabolism_rate = 2.5; % Rate at which Gluttons shrink
new_glutt_threshold=10;
%spawnRate = 0.1; % Rate at which Muttons spawn
muttonGirthIncrease = 5; % Girth increase per eaten Mutton
tsteps = 100; % Simulation time
spawn_rate= 10; %mutton spawn rate, currently unused
init_mutton = 50;
mutton_threshold=0;

% Initialize Gluttons
gluttons = repmat(struct('pos', []), 1, num_glutt);
for i = 1:num_glutt
    gluttons(i).pos = randi(grid_size, 1, 2); % Random position 
    gluttons(i).girth = initial_girth; % Initial size
end

% Init Muttons with flag
muttons = repmat(struct('pos', [], 'eaten', false), 1, init_mutton);
for j = 1:init_mutton
    muttons(j).pos = randi(grid_size, 1, 2); % Random position for each Mutton
    muttons(j).eaten = false; % Initialize as not eaten
end

% Init arrays to store population over time
glutton_Pop = zeros(1, tsteps);
mutton_Pop = zeros(1, tsteps);
% Set the init population counts
glutton_Pop(1) = num_glutt; % Initial number of Gluttons
mutton_Pop(1) = init_mutton; % Initial number of Muttons

% Simulation loop
for t = 2:tsteps
    % move Gluttons 
    for i = 1:num_glutt
        if i <= length(gluttons)
            move = randi([-1, 1], 1, 2); % Random move in x and y
            gluttons(i).pos = gluttons(i).pos + move;
            % gluttons within bounds
            gluttons(i).pos = max(min(gluttons(i).pos, grid_size), [1, 1]);
            % metabolism reduces size
            gluttons(i).girth = gluttons(i).girth - metabolism_rate;
        end    
        % Check if Glutton's girth is below 0 and remove 
        i = 1; % Start from the first glutton
        while i <= length(gluttons)
            if gluttons(i).girth <= 0.1
                gluttons(i) = []; % Remove Glutton
                num_glutt = num_glutt - 1; % Decrease the count of Gluttons
            else
                i = i + 1; % Move to the next Glutton
            end
        end    
    end   
    
    % Gluttons eat Muttons
    for i = 1:num_glutt
        j = 1; % Start with the first Mutton
        while j <= length(muttons)
            if i <= length(gluttons) && ~muttons(j).eaten %go through gluttons and non eaten muttons
                % Calculate distance between glutton and current Mutton
                distance = norm(gluttons(i).pos - muttons(j).pos);
                if distance <= 1 % distance threshold
                    % Mutton is eaten
                    disp(['time eaten:', num2str(t)]);
                    gluttons(i).girth = gluttons(i).girth + muttonGirthIncrease;
                    muttons(j).eaten = true; % Mark the Mutton as eaten
                    muttons(j) = []; % Remove eaten Mutton, including its position
                    disp(['Number of Muttons after removal: ', num2str(length(muttons))]);
                    % If a new Glutton should be spawned
                    if gluttons(i).girth >= new_glutt_threshold
                        num_glutt = num_glutt + 1; % Increase the count of Gluttons
                        gluttons(num_glutt).pos = gluttons(i).pos; % New Glutton starts at the same position
                        gluttons(num_glutt).girth = initial_girth; % New Glutton starts with initial girth
                        gluttons(i).girth = gluttons(i).girth - new_glutt_threshold; % Reset original Glutton's girth
                    end
                    break;
                else
                    j = j + 1; % Move to the next Mutton if the current one is not eaten
                end
            else
                break; % Break the loop if no Gluttons are left
            end
        end
    end

    % Determine whether to spawn Muttons based on time and current number of Muttons
    if length(muttons) > mutton_threshold
        % Compute spawn probability based on time step
        if mod(t, 1) == 0
            % Create a new Mutton with a random position and the 'eaten' field set to false
            newMutton = struct('pos', randi(grid_size, 1, 2), 'eaten', false);
            % Add the new Mutton to the array of Muttons
            muttons = [muttons, newMutton];
            
        end   
    end    
    %end
    % Update mutton population count at each time step
    glutton_Pop(t) = num_glutt; % gluttons over time
    mutton_Pop(t) = length(muttons); % The number of muttons
end %closes for t=2:tsteps   

% Plot the current state
clf;
% Plotting the population over time
plot(1:tsteps, glutton_Pop, 'r', 'LineWidth', 2);
hold on;
plot(1:tsteps, mutton_Pop, 'g', 'LineWidth', 2);
legend('glutton', 'mutton')
xlabel('Time');
ylabel('Population');
title('Population of Gluttons and Muttons Over Time');
grid on;
    