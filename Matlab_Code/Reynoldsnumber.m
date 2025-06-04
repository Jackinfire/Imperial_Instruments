clear all; close all; clc;

% Define range of Delta_P
Delta_P_values = 20000:10000:150000; % Pressure difference range (Pa)

% Channel widths to plot
channel_widths = [4e-5, 8e-5, 16e-5];
titles = {'Width = 4e-5 m', 'Width = 8e-5 m', 'Width = 16e-5 m'};
% Define custom colors for each plot (you can add more colors as needed)
custom_colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4940, 0.1840, 0.5560]}; 

% Create a figure for the line plot
fig = figure('Position', [100, 100, 1200, 900], 'Color', [1, 1, 1]);
hold on;

for k = 1:length(channel_widths)
    channel_width = channel_widths(k);

    % Initialize array to store maximum Reynolds numbers
    max_reynolds_numbers = zeros(size(Delta_P_values));

    % Loop through each Delta_P and calculate maximum Reynolds number
    for i = 1:length(Delta_P_values)
        Delta_P = Delta_P_values(i);
        max_reynolds_numbers(i) = computeMaxReynoldsNumber(Delta_P, channel_width); 
    end

    % Get the custom color for this plot (cycling through the list if needed)
    color_index = mod(k-1, length(custom_colors)) + 1;
    plot_color = custom_colors{color_index};

    % Plot the results on the same graph
    plot(Delta_P_values, max_reynolds_numbers, '-', 'LineWidth', 2, 'MarkerSize', 8, ...
        'Color', plot_color, 'DisplayName', titles{k});
end

% Enhanced plot features
xlabel('\textbf{Pressure Difference} $\Delta P$ \textbf{(Pa)}', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('\textbf{Reynolds Number}', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 14);
grid on;

% Customize axes
ax = gca;
set(ax, 'Color', [0.95, 0.95, 0.95]); % Light gray background for axes
set(ax, 'XColor', 'k', 'YColor', 'k'); % Black axis labels
set(ax, 'LineWidth', 1.5, 'FontSize', 14);
set(ax, 'GridColor', 'k', 'MinorGridColor', 'k', 'GridLineStyle', '--', 'MinorGridLineStyle', ':', 'MinorGridAlpha', 0.5);

% Adjust axis padding
outerpos = ax.OuterPosition;
ti = ax.TightInset;
padding = 0.05; % Adjust this value for more/less padding
ax.Position = [outerpos(1) + ti(1) + padding, outerpos(2) + ti(2) + padding, ...
               outerpos(3) - ti(1) - ti(3) - 2 * padding, ...
               outerpos(4) - ti(2) - ti(4) - 2 * padding];

% Store results in a struct array for saving
reynolds_data = struct();
reynolds_data.Delta_P_values = Delta_P_values;
reynolds_data.channel_widths = channel_widths;
reynolds_data.results = [];

for k = 1:length(channel_widths)
    reynolds_data.results(k).width = channel_widths(k);
    reynolds_data.results(k).title = titles{k};
    reynolds_data.results(k).color = custom_colors{k};
    reynolds_data.results(k).max_reynolds_numbers = zeros(size(Delta_P_values));

    for i = 1:length(Delta_P_values)
        Delta_P = Delta_P_values(i);
        reynolds_data.results(k).max_reynolds_numbers(i) = computeMaxReynoldsNumber(Delta_P, channel_widths(k));
    end
end

hold off;
% Save the computed Reynolds data to a .mat file
save('reynolds_vs_pressure.mat', 'reynolds_data');


% Export the figure with high resolution and proper size
% exportgraphics(gcf, '/Users/ommahajan/Desktop/Year_3/BIOE60005 - Bioengineering Group Project/Report Stuff/Reynolds.png', ...
  %  'BackgroundColor', 'white', 'Resolution', 600);

% Function to compute maximum Reynolds number
function max_reynolds = computeMaxReynoldsNumber(Delta_P, channel_width)
    % Parameters
    channel_height = 4e-5;   % Channel height (m)
    mu = 0.001;              % Dynamic viscosity (Pa·s)
    rho = 1000;              % Density of water (kg/m^3)
    L = 1e-2;                % Characteristic length (m)
    N_terms = 50;            % Number of terms in the summation

    % Derived parameters
    h = channel_height;    % Full channel height (m)
    w = channel_width;     % Channel width (m)

    % Define spatial ranges
    y_range = linspace(-w/2, w/2, 100); % y-coordinates (m)
    z_range = linspace(0, h, 100);      % z-coordinates (m)

    % Initialize velocity matrix
    u = zeros(length(y_range), length(z_range));

    % Compute velocity profile
    for i = 1:length(y_range)
        for j = 1:length(z_range)
            y = y_range(i);
            z = z_range(j);
            % Summation over odd n
            sum_term = 0;
            for n = 1:2:(2 * N_terms - 1)
                cosh_term = cosh((n * pi * y) / h) / cosh((n * pi * w) / (2 * h));
                sin_term = sin((n * pi * z) / h);
                sum_term = sum_term + (1 / n^3) * (1 - cosh_term) * sin_term;
            end
            % Velocity at (y, z)
            u(i, j) = ((4 * h^2 * Delta_P) / (pi^3 * mu * L)) * sum_term; % Velocity in m/s
        end
    end

    % Extract velocity at z = h/2 (j = 50)
    j = 50; 
    velocity_at_mid_z = u(:, j);
    
    % Maximum velocity
    max_velocity = max(abs(velocity_at_mid_z));
    
    % Hydraulic diameter
    D_h = (2 * w * h) / (w + h);
    
    % Maximum Reynolds number calculation
    max_reynolds = (rho * max_velocity * D_h) / mu;

end
