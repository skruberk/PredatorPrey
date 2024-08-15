%mutton and glutton, predator-prey

%when gluttons eat a mutton their size increases by a variable amt., when size
%reaches threshold you spawn a new mutton
%muttons spawn randomly as long as they're above a threshold population
%muttons stay put while gluttons move
%glutton movement imparts penalty to their size (proxy for metabolism)

clc;    % Clear the command window.
clearvars;
close all;  % Close all figs
workspace;  % Make sure the workspace panel is showing.
% Parameters
num_glutt = 5; % Number of gluttons
grid_size = [100, 100]; % size of the 2D world
initial_girth = 5; % initial size of gluttons
metabolism_rate = 1.5; % glutton shrink rate, movement penalty
new_glutt_threshold=5;
%spawnRate = 0.1; % Rate at which Muttons spawn, now runs a different way
muttongirthinc = 2.5; % size increase of glutton per eaten Mutton
tsteps = 50; % Simulation time
%spawn_rate= 10; %mutton spawn rate, currently unused, now calculated based on time
init_mutton = 50; %initial mutton population 
mutton_threshold=0; % check for positive mutton values 

% initialize Gluttons
gluttons = repmat(struct('pos', []), 1, num_glutt);
for i = 1:num_glutt
    gluttons(i).pos = randi(grid_size, 1, 2); % Random position 
    gluttons(i).girth = initial_girth; % Initial size
end

% init muttons with flag
muttons = repmat(struct('pos', [], 'eaten', false), 1, init_mutton);
for j = 1:init_mutton
    muttons(j).pos = randi(grid_size, 1, 2); % Random position for each Mutton
    muttons(j).eaten = false; % Initialize as not eaten
end

% init arrays to store population over time
glutton_Pop = zeros(1, tsteps);
mutton_Pop = zeros(1, tsteps);
% Set the init population counts
glutton_Pop(1) = num_glutt; % Initial number of Gluttons
mutton_Pop(1) = init_mutton; % Initial number of Muttons

% simulation loop
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
                gluttons(i) = []; % remove glutton
                num_glutt = num_glutt - 1; % Decrease the count of gluttons
            else
                i = i + 1; % move to the next glutton
            end
        end    
    end   
    
    % gluttons eat muttons
    for i = 1:num_glutt
        j = 1; % Start with the first Mutton
        while j <= length(muttons)
            if i <= length(gluttons) && ~muttons(j).eaten %go through gluttons and non eaten muttons
                % Calculate distance between glutton and current Mutton
                distance = norm(gluttons(i).pos - muttons(j).pos);
                if distance <= 1 % distance threshold
                    % Mutton is eaten
                    disp(['time eaten:', num2str(t)]);
                    gluttons(i).girth = gluttons(i).girth + muttongirthinc;
                    muttons(j).eaten = true; % mark the Mutton as eaten
                    muttons(j) = []; % Remove eaten Mutton, including its position
                    disp(['Number of Muttons after removal: ', num2str(length(muttons))]);
                    % if a new Glutton should be spawned..
                    if gluttons(i).girth >= new_glutt_threshold
                        num_glutt = num_glutt + 1; % increase the count of Gluttons
                        gluttons(num_glutt).pos = gluttons(i).pos; % new glutton starts at the same position
                        gluttons(num_glutt).girth = initial_girth; % new glutton starts with initial girth
                        gluttons(i).girth = gluttons(i).girth - new_glutt_threshold; % reset to original glutton girth
                    end
                    break;
                else
                    j = j + 1; % move to next Mutton if the current one is not eaten
                end
            else
                break; % Break the loop if no Gluttons are left
            end
        end
    end %closes for i:num_glutt

    % spawn muttons based on time and current number of Muttons
    if length(muttons) > mutton_threshold
        % Compute spawn probability based on time step
        if mod(t, 1) == 0
            % Randomly decide whether to spawn 1 or 2 Muttons
            numNewMuttons = randi([1, 10]);
            for i = 1:numNewMuttons
                % Create a new Mutton with a random position and the 'eaten' field set to false
                newMutton = struct('pos', randi(grid_size, 1, 2), 'eaten', false);
                % Add the new Mutton to the array of Muttons
                muttons = [muttons, newMutton];
            end
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
title('Population of Gluttons & Muttons Over Time');
grid on;
    