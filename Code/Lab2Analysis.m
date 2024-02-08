% AER E 344 Spring 2024 Lab 02 Analysis
% Section 3 Group 3
clear, clc, close all;

figure_dir = "../Figures/";
u = symunit;

%% Import Data
data_sheet = readtable('AER E 344 Lab 02 Data Sheet.xlsx', ...
    'VariableNamingRule', 'preserve');
omega_motor = data_sheet.("Motor speed [Hz]").'; % [Hz]
H_ref = double(separateUnits(unitConvert( ...
    data_sheet.("H_ref [in.]").' * u.in, u.m))); % [m]
H_A = double(separateUnits(unitConvert( ...
    data_sheet.("H_A [in.]").' * u.in, u.m))); % [m]
H_E = double(separateUnits(unitConvert( ...
    data_sheet.("H_E [in.]").' * u.in, u.m))); % [m]
H_total = double(separateUnits(unitConvert( ...
    data_sheet.("H_total [in.]").' * u.in, u.m))); % [m]
H_static = double(separateUnits(unitConvert( ...
    data_sheet.("H_static [in.]").' * u.in, u.m))); % [m]
T_tunnel = data_sheet.("T_tunnel [deg C]").'; % [ºC]

%% Variables
% https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html
% Calculated @ 22.1ºC
rho_water = 997.74; % [kg / m^3]
% https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html
% Calculated @ 22.1ºC
rho_air = 1.195; % [kg / m^3]
% https://physics.nist.gov/cgi-bin/cuu/Value?gn
g = 9.80665; % [m / s^2]

%% Calculate q_T & delta_p
% q_T = P_0T - P_T
q_T = rho_water .* g .* (H_static - H_total); % [Pa]
delta_p = rho_water .* g .* (H_E - H_A); % [Pa]

%% Calculate q_T vs. delta_p Regression
[regress_1, S_1] = polyfit(delta_p, q_T, 1);
K = regress_1(1); % []
regress_1_x = delta_p(1):0.1:delta_p(end); % [Pa]
regress_1_y = polyval(regress_1, regress_1_x); % [Pa]

SStot_1 = sum((q_T - mean(q_T)).^2);
SSres_1 = sum((q_T - polyval(regress_1, delta_p)).^2);
Rsq_1 = 1 - SSres_1 / SStot_1; % R^2 value

fprintf("K = %g []\n" + ...
    "R^2 (q_T vs. delta_p) = %g\n", ...
    K, Rsq_1);

%% Plot q_t vs delta_p
figure(1);
scatter(delta_p, q_T);
title("Dynamic Pressure vs Change in Static Pressure")
xlabel("{\Delta}p [Pa]")
ylabel("q_t [Pa]")
hold on;
plot(regress_1_x, regress_1_y);
hold off;
legend("Experimental Data", "Line of Best Fit", "Location", "northwest");
grid on;
saveas(gcf, ...
    figure_dir + "Dynamic Pressure vs Change in Static Pressure.svg");

%% Calculate v_T
v_T = sqrt(2 * q_T / rho_air); % [m/s]

%% Calculate v_T vs omega_motor Regression
[regress_2, S_2] = polyfit(omega_motor, v_T, 1);
regress_2_x = omega_motor(1):0.1:omega_motor(end); % [Hz]
regress_2_y = regress_2(1) * regress_2_x + regress_2(2); % [m/s]

SStot_2 = sum((v_T - mean(v_T)).^2);
SSres_2 = sum((v_T - polyval(regress_2, omega_motor)).^2);
Rsq_2 = 1 - SSres_2 / SStot_2; % R^2 value

fprintf("v_T = %G * omega_motor + %g\n" + ...
    "R^2 (v_T vs omega_motor) = %g\n", ...
    regress_2, Rsq_2);

%% Plot v_T vs omega_motor
figure(2);
scatter(omega_motor, v_T);
title("Test Chamber Velocity vs Motor Frequency");
xlabel("\omega_{motor} [Hz]");
ylabel("v_T [m/s]");
hold on;
plot(regress_2_x, regress_2_y);
hold off;
legend("Experimental Data", "Line of Best Fit", "Location", "northwest");
grid on;
saveas(gcf, figure_dir + "Test Chamber Velocity vs Motor Frequency.svg");
