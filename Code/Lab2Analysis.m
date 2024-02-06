clear, clc, close all;

u = symunit;

%% Import Data
Data_Sheet = readtable('AER E 344 Lab 02 Data Sheet.xlsx', ...
    'VariableNamingRule', 'preserve');
Motor_fr = Data_Sheet.("Motor speed [Hz]").'; % [Hz]
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

%% Calculate q_T & delta_P
% q_T = P_0T - P_T
q_T = rho_water * g * (H_static - H_total); % [Pa]
delta_P = rho_water * g * (H_E - H_A); % [Pa]

K_array = q_T ./ delta_P;

v_T = sqrt((2 .* K_array .*delta_P)/rho_air);
q = 0.5 * rho_air .* v_T.^2; 

%% Plot
f1 = figure('Name', 'Motor Freq vs Air Velocity');
hold on
v_T(1) = 0;
scatter(Motor_fr, v_T);
xlabel("Motor Frequency (Hz)")
ylabel("Air Velocity (m/s)")
regress = polyfit(Motor_fr, v_T, 1);
plot(Motor_fr, Motor_fr*regress(1));
fprintf("Motor Freq vs Air Vel Slope = %f\n", regress(1))

f2 = figure('Name', 'Delta P vs q');
hold on
q(1)=0;
scatter(delta_P, q);
xlabel("Delta P")
ylabel(".5*p*V^2")
K = polyfit(delta_P,q,1);
plot(delta_P, delta_P*K(1));
fprintf("K = %f\n", K(1))
