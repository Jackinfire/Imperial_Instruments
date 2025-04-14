clear all; close all;

% Define range of Delta_P
Delta_P_values = 20000:2000:150000; % Pressure difference range (Pa)

% Channel widths to plot
channel_widths = [4e-5, 8e-5, 16e-5];
titles = {'Width = 4e-5 m', 'Width = 8e-5 m', 'Width = 16e-5 m'};
custom_colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4940, 0.1840, 0.5560]}; 

% Create a single figure for combined plot
fig = figure('Position', [100, 100, 1200, 900], 'Color', [1, 1, 1]);
hold on;

% Initialize struct to store all data
shear_data = struct();
shear_data.Delta_P_values = Delta_P_values;
shear_data.channel_widths = channel_widths;
shear_data.results = [];

for k = 1:length(channel_widths)
    channel_width = channel_widths(k);
    shear_data.results(k).width = channel_width;
    shear_data.results(k).title = titles{k};
    shear_data.results(k).color = custom_colors{k};
    
    % Initialize array to store maximum shear rates
    max_shear_rates = zeros(size(Delta_P_values));
    
    for i = 1:length(Delta_P_values)
        Delta_P = Delta_P_values(i);
        max_shear_rates(i) = computeMaxShearRate(Delta_P, channel_width); 
    end
    
    shear_data.results(k).max_shear_rates = max_shear_rates;
    
    % Plot
    color_index = mod(k-1, length(custom_colors)) + 1;
    plot_color = custom_colors{color_index};
    plot(Delta_P_values, max_shear_rates, '-.', 'LineWidth', 2, 'MarkerSize', 5, ...
        'Color', plot_color, 'DisplayName', titles{k});
end

% Enhanced plot features
xlabel('\textbf{Pressure Difference} $\Delta P$ \textbf{(Pa)}', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('\textbf{Maximum Shear Rate (1/s)}', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
% Customize axes
ax = gca;
set(ax, 'Color', [0.95, 0.95, 0.95]); % Light gray background for axes
set(ax, 'XColor', 'k', 'YColor', 'k'); % Black axis labels
set(ax, 'LineWidth', 1.5, 'FontSize', 14);
grid on;
legend('Location', 'best', 'FontSize', 14);

hold off;
save('shear_rate_vs_pressure.mat', 'shear_data');


% Function to compute maximum shear rate
function max_gradient = computeMaxShearRate(Delta_P, channel_width)
    % Parameters
    channel_height = 4e-5;   % Channel height (m)
    mu = 0.001;              % Dynamic viscosity (Pa·s)
    L = 1e-2;                % Characteristic length (m)
    N_terms = 50;            % Number of terms in the summation

    % Derived parameters
    h = channel_height;    % Full channel height (m)
    w = channel_width;     % Channel width (m)

    % Define spatial ranges
    y_range = linspace(-w/2, w/2, 100); % y-coordinates (m)
    z_range = linspace(0, h, 100);      % z-coordinates (m)
    dy = y_range(2) - y_range(1);       % Step size in y

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
    
    % Gradient with respect to y
    gradient_velocity = diff(velocity_at_mid_z) / dy;
    
    % Maximum gradient
    max_gradient = max(abs(gradient_velocity));
end
