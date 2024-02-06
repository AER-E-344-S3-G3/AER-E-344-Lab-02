% AER E 344 Spring 2024 Lab 02 Analysis
% Section 3 Group 3
clear, clc, close all;

u = symunit;

%% Import Data
Data_Sheet = readtable('AER E 344 Lab 02 Data Sheet.xlsx', ...
    'VariableNamingRule', 'preserve');
omega_motor = Data_Sheet.("Motor speed [Hz]").'; % [Hz]
H_ref = double(separateUnits(unitConvert( ...
    Data_Sheet.("H_ref [in.]").' * u.in, u.m))); % [m]
H_A = double(separateUnits(unitConvert( ...
    Data_Sheet.("H_A [in.]").' * u.in, u.m))); % [m]
H_E = double(separateUnits(unitConvert( ...
    Data_Sheet.("H_E [in.]").' * u.in, u.m))); % [m]
H_total = double(separateUnits(unitConvert( ...
    Data_Sheet.("H_total [in.]").' * u.in, u.m))); % [m]
H_static = double(separateUnits(unitConvert( ...
    Data_Sheet.("H_static [in.]").' * u.in, u.m))); % [m]
T_tunnel = Data_Sheet.("T_tunnel [deg C]").'; % [ºC]

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

%% Calculate q_t vs. delta_p Regression
regress_1 = polyfit(delta_p, q_T, 1);
K = regress_1(1); % []
regress_1_x = delta_p(1):0.1:delta_p(end); % [Pa]
regress_1_y = K * regress_1_x + regress_1(2); % [Pa]

fprintf("K = %g []\n", K);

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

%% Calculate v_T
v_T = sqrt(2 * q_T / rho_air); % [m/s]

%% Calculate v_T vs omega_motor Regression
regress_2 = polyfit(omega_motor, v_T, 1);
regress_2_x = omega_motor(1):0.1:omega_motor(end); % [Hz]
regress_2_y = regress_2(1) * regress_2_x + regress_2(2); % [m/s]

fprintf("v_T = %G * omega_motor + %g\n", regress_2);

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
