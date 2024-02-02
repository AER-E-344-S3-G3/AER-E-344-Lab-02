clear, clc, close all

in_2_m = 0.0254;

%% Data Sheet
Data_Sheet = readtable('AER E 344 Lab 02 Data Sheet.xlsx','VariableNamingRule','preserve');
Motor_fr = Data_Sheet.("Motor speed [Hz]").';
H_ref = Data_Sheet.("H_ref [in.]").' * in_2_m;
H_A = Data_Sheet.("H_A [in.]").' * in_2_m;
H_E = Data_Sheet.("H_E [in.]").' * in_2_m;
H_total = Data_Sheet.("H_total [in.]").' * in_2_m;
H_static = Data_Sheet.("H_static [in.]").' * in_2_m;
T_tunnel = Data_Sheet.("T_tunnel [deg C]").' * in_2_m;

%% Variables
P_ref = 100605; % N/m^2 .9929atm * 101325
p_water = 997.77; % kg / m^3
p_air = 1.225; % kg/m^3
g = 9.8; % m/s^2

%% Calculate q_T & delta_p at each data point
P_01T = P_ref - p_water * g * (H_total - H_ref);
P_T = P_ref - p_water * g * (H_static - H_ref);
q_T = P_01T - P_T;

P_A = P_ref - p_water * g * (H_A - H_ref);
P_E = P_ref - p_water * g * (H_E - H_ref);
delta_p = P_A - P_E;

K_array = q_T ./ delta_p;

v_T = sqrt((2 .* K_array .*delta_p)/p_air);
q = 0.5 * p_air .* v_T.^2; 

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
scatter(delta_p, q);
xlabel("Delta P")
ylabel(".5*p*V^2")
K = polyfit(delta_p,q,1);
plot(delta_p, delta_p*K(1));
fprintf("K = %f\n", K(1))

