% Script to evaluate Simulation2 for different N_terms
% and determine optimal N_terms based on accuracy and computation time

% Input parameters (SI units)
channel_height = 4e-5; % Channel height (m)
channel_width = 16e-5;  % Channel width (m)
Delta_P = 10000;         % Pressure difference (Pa)
mu = 0.001;            % Dynamic viscosity (Pa·s)
L = 1e-2;                 % Characteristic length (m)

% Range of N_terms to test
N_values = 15:5:100;
execution_times = zeros(size(N_values));
differences = zeros(size(N_values) - [0, 1]); % Store differences between profiles

% Reference grid for velocity profile comparison
y_range = linspace(-channel_width/2, channel_width/2, 100); % y-coordinates (m)
z_range = linspace(0, channel_height, 100);                % z-coordinates (m)
[Y, Z] = meshgrid(y_range, z_range);

% Loop through N_values
prev_profile = []; % Placeholder for previous velocity profile
for idx = 1:length(N_values)
    N_terms = N_values(idx);
    
    fprintf('Running Simulation2 with N_terms = %d...\n', N_terms);
    
    % Measure computation time
    tic;
    % Call Simulation2 but only compute the profile (no plotting)
    u = compute_velocity_profile(channel_height, channel_width, Delta_P, mu, L, N_terms, Y, Z);
    execution_times(idx) = toc;
    
    % Calculate difference if not the first run
    if ~isempty(prev_profile)
        diff = abs(u - prev_profile); % Absolute difference
        differences(idx - 1) = max(diff(:)); % Store max difference as accuracy metric
    end
    
    % Update previous profile
    prev_profile = u;
    
end
% Create figure with specified size and white background color
fig = figure('Position', [100, 100, 1200, 900], 'Color', [1, 1, 1]); 
hold on;

% Improved line plot
plot(N_values(2:end), differences, '-x', 'LineWidth', 2, 'MarkerSize', 10, ...
    'Color', [0, 0.4470, 0.7410], 'MarkerEdgeColor', [0.8, 0.2, 0.2]);

% Enhanced plot features with bold axes titles
xlabel('\textbf{Number of Terms}', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
ylabel('\textbf{Max Difference}', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
grid on;

% Adjusting axis limits and grid appearance
xlim([min(N_values) max(N_values)]);
set(gca, 'GridLineStyle', '--', 'LineWidth', 1.5, 'FontSize', 16, ...
    'Color', [0.9, 0.9, 0.9], ...       % Light gray background for the inner box
    'XColor', 'k', 'YColor', 'k', ...   % Black axis lines and labels
    'GridColor', 'k', 'MinorGridColor', 'k'); % Black grid lines

% Minor grid lines for better visibility
set(gca, 'MinorGridLineStyle', ':', 'MinorGridAlpha', 0.5);

% Adjusting axis position for padding
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
% Add padding by increasing the outer position margins
padding = 0.05; % Adjust this value for more/less padding
ax.Position = [outerpos(1) + ti(1) + padding, outerpos(2) + ti(2) + padding, ...
               outerpos(3) - ti(1) - ti(3) - 2 * padding, ...
               outerpos(4) - ti(2) - ti(4) - 2 * padding];

hold off;

% Nested function to compute velocity profile without plotting
function u = compute_velocity_profile(h, w, Delta_P, mu, L, N_terms, Y, Z)
    % Initialize velocity matrix
    u = zeros(size(Y));
    
    % Compute velocity profile
    for i = 1:size(Y, 1)
        for j = 1:size(Y, 2)
            y = Y(i, j);
            z = Z(i, j);
            % Summation over odd n
            sum_term = 0;
            for n = 1:2:(2*N_terms-1)
                cosh_term = cosh((n*pi*y)/h) / cosh((n*pi*w)/(2*h));
                sin_term = sin((n*pi*z)/h);
                sum_term = sum_term + (1/n^3) * (1 - cosh_term) * sin_term;
            end
            % Velocity at (y, z)
            u(i, j) = (4*h^2*Delta_P) / (pi^3*mu*L) * sum_term; % Velocity in m/s
        end
    end
end
