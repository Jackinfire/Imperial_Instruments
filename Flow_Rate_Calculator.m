clear all;close all;clc;

% Define constants
h = 4e-5; % Height of the channel (in meters)
w = 12e-5;  % Width of the channel (in meters)
L = 0.01;  % Length of the channel (in meters)
mu = 0.001; % Dynamic viscosity (in Pascal-seconds)
deltaP = linspace(9000,140000, 1000); % Pressure drop (in Pascals), variable from 0 to 10 Pa

% Calculate Q for each deltaP
Q = zeros(size(deltaP)); % Initialize flow rate array
N = 60; % Number of terms in the series expansion

for i = 1:length(deltaP)
    sum_term = 0; % Initialize summation term
    for n = 1:2:2*N % Only odd values of n
        sum_term = sum_term + (1/n^5) * (192 * h / (pi^5 * w)) * tanh(n * pi * w / (2 * h));
    end
    Q(i) = (h^3 * w * deltaP(i)) / (12 * mu * L) * (1 - sum_term);
end

% Convert Q from cubic meters per second to milliliters per hour
Q_ml_per_hr = Q * 1e6 * 3600;

% Plot the results
figure;
plot(deltaP, Q_ml_per_hr, 'LineWidth', 2);
xlabel('\Delta P (Pa)', 'FontSize', 12);
ylabel('Q (mL/hr)', 'FontSize', 12);
title('Volumetric Flow Rate vs Pressure Drop', 'FontSize', 14);
grid on;