function [name, PID, dataOut1, dataOut2] = Wing_Analysis_Function(dataIn)
dataIn = [201,2.1500E+02,6.0000E+01,9.6000E+00,5.0000E+00,5.0000E-02,1.0000E-01,3.7500E+00,2.4000E+01,2.8000E+01,0.0000E+00,6.0000E-01,5.2100E-03,5.2100E-03,0.0000E+00,1.0000E-01,1.0000E+01,3.7000E+01,4.3000E+01,-3.7000E+01,-4.3000E+01,3.0000E+01,4.2000E+00,1.0400E-02,4.1700E-02,0.0000E+00,1.0000E-01,1.0000E+01,3.7000E+01,4.3000E+01,-3.7000E+01,-4.3000E+01,6.0000E+01,1.5000E+00,5.2100E-03,5.2100E-03,0.0000E+00,1.0000E-01,1.0000E+01,3.7000E+01,4.3000E+01,-3.7000E+01,-4.3000E+01,1.5000E+01,3.5000E+00,1.0400E-02,4.1700E-02,0.0000E+00,1.0000E-01,1.0000E+01,3.7000E+01,4.3000E+01,-3.7000E+01,-4.3000E+01,1.1000E+00,1.5000E+00,3.8000E+00,2.4000E+00,4.0000E-01,1.0000E+01,5.7000E+01,-1.5200E+01,-3.8000E+00,1.5000E+01,0.0000E+00];

% Analysis Parameters
nplot    = 201;   % number of output plot data points
Lo       = 2.1500E+02;   % Wing Length (inch)
Co       = 6.0000E+01;   % Wing Chord (inch)
tmax     = 9.6000E+00;   % Maximum Wing Thickness (inch)
Kaw      = 5.0000E+00;   % Secondary Structure Added Weight (%)
to_sk    = 5.0000E-02;   % Wing Skin Thickness (inch)
rho_sk   = 1.0000E-01;   % Wing Skin Weight Density (lb/inch^3)
Go_sk    = 3.7500E+00;   % Skin Material Shear Modulus (Msi)
Sys_sk   = 2.4000E+01;   % Skin Material Yield Shear Strength (Ksi)
Sus_sk   = 2.8000E+01;   % Skin Material Ultimate Shear Strength (Ksi)
yo_str1  = 0.0000E+00;   % Stringer 1 y-location (inch)
A_str1   = 6.0000E-01;   % Stringer 1 Cross-Section Area (inch^2)
Iyy_str1 = 5.2100E-03;   % Stringer 1 Iyy Inertia about the y-axis (inch^4)
Izz_str1 = 5.2100E-03;   % Stringer 1 Izz Inertia about the z-axis (inch^4)
Iyz_str1 = 0.0000E+00;   % Stringer 1 Iyz Product of Inertia (inch^4)
rho_str1 = 1.0000E-01;   % Stringer 1 Wing Skin Weight Density (lb/inch^3)
Eo_str1  = 1.0000E+01;   % Stringer 1 Material Young's Modulus (Msi)
Syt_str1 = 3.7000E+01;   % Stringer 1 Material Yield Strength    - Tension (Ksi)
Sut_str1 = 4.3000E+01;   % Stringer 1 Material Ultimate Strength - Tension (Ksi)
Syc_str1 = -3.7000E+01;   % Stringer 1 Material Yield Strength    - Compression (Ksi)
Suc_str1 = -4.3000E+01;   % Stringer 1 Material Ultimate Strength - Compression (Ksi)
yo_str2  = 3.0000E+01;   % Stringer 2 y-location (inch)
A_str2   = 4.2000E+00;   % Stringer 2 Cross-Section Area (inch^2)
Iyy_str2 = 1.0400E-02;   % Stringer 2 Iyy Inertia about the y-axis (inch^4)
Izz_str2 = 4.1700E-02;   % Stringer 2 Izz Inertia about the z-axis (inch^4)
Iyz_str2 = 0.0000E+00;   % Stringer 2 Iyz Product of Inertia (inch^4)
rho_str2 = 1.0000E-01;   % Stringer 2 Wing Skin Weight Density (lb/inch^3)
Eo_str2  = 1.0000E+01;   % Stringer 2 Material Young's Modulus (Msi)
Syt_str2 = 3.7000E+01;   % Stringer 2 Material Yield Strength    - Tension (Ksi)
Sut_str2 = 4.3000E+01;   % Stringer 2 Material Ultimate Strength - Tension (Ksi)
Syc_str2 = -3.7000E+01;   % Stringer 2 Material Yield Strength    - Compression (Ksi)
Suc_str2 = -4.3000E+01;   % Stringer 2 Material Ultimate Strength - Compression (Ksi)
yo_str3  = 6.0000E+01;   % Stringer 3 y-location (inch)
A_str3   = 1.5000E+00;   % Stringer 3 Cross-Section Area (inch^2)
Iyy_str3 = 5.2100E-03;   % Stringer 3 Iyy Inertia about the y-axis (inch^4)
Izz_str3 = 5.2100E-03;   % Stringer 3 Izz Inertia about the z-axis (inch^4)
Iyz_str3 = 0.0000E+00;   % Stringer 3 Iyz Product of Inertia (inch^4)
rho_str3 = 1.0000E-01;   % Stringer 3 Wing Skin Weight Density (lb/inch^3)
Eo_str3  = 1.0000E+01;   % Stringer 3 Material Young's Modulus (Msi)
Syt_str3 = 3.7000E+01;   % Stringer 3 Material Yield Strength    - Tension (Ksi)
Sut_str3 = 4.3000E+01;   % Stringer 3 Material Ultimate Strength - Tension (Ksi)
Syc_str3 = -3.7000E+01;   % Stringer 3 Material Yield Strength    - Compression (Ksi)
Suc_str3 = -4.3000E+01;   % Stringer 3 Material Ultimate Strength - Compression (Ksi)
yo_str4  = 1.5000E+01;   % Stringer 4 y-location (inch)
A_str4   = 3.5000E+00;   % Stringer 4 Cross-Section Area (inch^2)
Iyy_str4 = 1.0400E-02;   % Stringer 4 Iyy Inertia about the y-axis (inch^4)
Izz_str4 = 4.1700E-02;   % Stringer 4 Izz Inertia about the z-axis (inch^4)
Iyz_str4 = 0.0000E+00;   % Stringer 4 Iyz Product of Inertia (inch^4)
rho_str4 = 1.0000E-01;   % Stringer 4 Wing Skin Weight Density (lb/inch^3)
Eo_str4  = 1.0000E+01;   % Stringer 4 Material Young's Modulus (Msi)
Syt_str4 = 3.7000E+01;   % Stringer 4 Material Yield Strength    - Tension (Ksi)
Sut_str4 = 4.3000E+01;   % Stringer 4 Material Ultimate Strength - Tension (Ksi)
Syc_str4 = -3.7000E+01;   % Stringer 4 Material Yield Strength    - Compression (Ksi)
Suc_str4 = -4.3000E+01;   % Stringer 4 Material Ultimate Strength - Compression (Ksi)
SFy      = 1.1000E+00 ;   % Safety Factor - Yield
SFu      = 1.5000E+00 ;   % Safety Factor - Ultimate
LF       = 3.8000E+00 ;   % Aircraft Load Factor
py0      = 2.4000E+00;   % Drag Distribution - Constant (lb/in)
pyr      = 4.0000E-01;   % Drag Distribution - rth order (lb/in)
rth      = 1.0000E+01;   % Drag Distribution - polynomial order
pz0      = 5.7000E+01;   % Lift Distribution - Constant (lb/in)
pz2      = -1.5200E+01;   % Lift Distribution - 2nd Order (lb/in)
pz4      = -3.8000E+00;   % Lift Distribution - 4th Order (lb/in)
mx0      = 1.5000E+01;   % Twist Moment Distribution - Constant (lb-in/in)
mx1      = 0.0000E+00;   % Twist Moment Distribution - 1st Order (lb-in/in)

% Total Half-Span Wing Weight (W)
W_structure = 85.96;  % Weight of structural components, you need to calculate this based on your data
W_skin = rho_sk * Co * to_sk * Lo;  % Weight of wing skin
W_stringers = (rho_str1 * A_str1 + rho_str2 * A_str2 + rho_str3 * A_str3 + rho_str4 * A_str4) * Lo;  % Weight of stringers
W = W_structure + W_skin + W_stringers;

% Modulus Weighted Centroid (yc and zc)
yc = (yo_str1 * Eo_str1 * A_str1 + yo_str2 * Eo_str2 * A_str2 + yo_str3 * Eo_str3 * A_str3 + yo_str4 * Eo_str4 * A_str4) / (Eo_str1 * A_str1 + Eo_str2 * A_str2 + Eo_str3 * A_str3 + Eo_str4 * A_str4);
zc = ((yo_str1 * Eo_str1 * A_str1 + yo_str2 * Eo_str2 * A_str2 + yo_str3 * Eo_str3 * A_str3 + yo_str4 * Eo_str4 * A_str4) / (Eo_str1 * A_str1 + Eo_str2 * A_str2 + Eo_str3 * A_str3 + Eo_str4 * A_str4)) - 27.05514;

% Cross-Section Weight (rho_A)
rho_A = rho_str1 * A_str1 + 0.6998 + rho_str2 * A_str2 + rho_str3 * A_str3 + rho_str4 * A_str4;

% Axial Stiffness (EA)
EA = Eo_str1 * A_str1 + 97999902 + Eo_str2 * A_str2 + Eo_str3 * A_str3 + Eo_str4 * A_str4;

% Bending Stiffness (EIyy, EIzz, EIyz)
EIyy = Eo_str1 * Iyy_str1 + 1762899999.6878 + Eo_str2 * Iyy_str2 + Eo_str3 * Iyy_str3 + Eo_str4 * Iyy_str4;
EIzz = Eo_str1 * Izz_str1 + Eo_str2 * Izz_str2 + 26111999999.062 + Eo_str3 * Izz_str3 + Eo_str4 * Izz_str4;
EIyz = Eo_str1 * Iyz_str1 - 2607400000 + Eo_str2 * Iyz_str2 + Eo_str3 * Iyz_str3 + Eo_str4 * Iyz_str4;

% Torsional Stiffness (GJ)
GJ_Mu = 243334999.97073;
GJ_str1 = Go_sk * Iyy_str1;
GJ_str2 = Go_sk * Iyy_str2;
GJ_str3 = Go_sk * Iyy_str3;
GJ_str4 = Go_sk * Iyy_str4;

GJ = GJ_str1 + GJ_str2 + GJ_str3 + GJ_str4 + 4 * GJ_Mu;

% Shear Centers (ey, ez)
ey = ((yo_str1 * Go_sk * A_str1 + yo_str2 * Go_sk * A_str2 + yo_str3 * Go_sk * A_str3 + yo_str4 * Go_sk * A_str4) / (Go_sk * A_str1 + Go_sk * A_str2 + Go_sk * A_str3 + Go_sk * A_str4)) - 10.639;
ez = ((yo_str1 * Go_sk * A_str1 + yo_str2 * Go_sk * A_str2 + yo_str3 * Go_sk * A_str3 + yo_str4 * Go_sk * A_str4) / (Go_sk * A_str1 + Go_sk * A_str2 + Go_sk * A_str3 + Go_sk * A_str4)) - 27.5868;

% Display the results
disp(['Total Half-Span Wing Weight (W): ', num2str(W)]);
disp(['Modulus Weighted Centroid (yc, zc): ', num2str(yc), ', ', num2str(zc)]);
disp(['Cross-Section Weight (rho_A): ', num2str(rho_A)]);
disp(['Axial Stiffness (EA): ', num2str(EA)]);
disp(['Bending Stiffness (EIyy, EIzz, EIyz): ', num2str(EIyy), ', ', num2str(EIzz), ', ', num2str(EIyz)]);
disp(['Torsion Stiffness (GJ): ', num2str(GJ)]);
disp(['Shear Centers (ey, ez): ', num2str(ey), ', ', num2str(ez)]);

% Calculate Wing Cross-Section Perimeter (S)
S = Co + 2 * sqrt((Lo/2)^2 + (Co/2)^2);

% Calculate Additional wing weight from secondary or tertiary structures
Kw = Kaw / 100;

% Calculate Weight (nw)
nw = (S * rho_sk * to_sk + (rho_str1 * A_str1 + rho_str2 * A_str2 + rho_str3 * A_str3 + rho_str4 * A_str4) * Lo * (1 + Kw));

% Graph for Distributed Drag (lb/inch)
x = linspace(0, Lo/2, nplot);
py = py0 + pyr * (x / Lo).^rth;
figure;
plot(x, py, 'b', 'LineWidth', 2);
xlabel('Distance along wing span (inch)');
ylabel('Distributed Drag (lb/inch)');
title('Distributed Drag along the Wing Span');

% Graph for Distributed Lift (lb/inch)
pz = pz0 + pz2 * (x/Lo).^2 + pz4 * (x/Lo).^4;
figure;
plot(x, pz, 'r', 'LineWidth', 2);
xlabel('Distance along wing span (inch)');
ylabel('Distributed Lift (lb/inch)');
title('Distributed Lift along the Wing Span');

% Graph for Distributed Aerodynamic Torque (lb-inch/inch)
Mx = mx0 + mx1 * (x/Lo);
figure;
plot(x, Mx, 'g', 'LineWidth', 2);
xlabel('Distance along wing span (inch)');
ylabel('Distributed Aerodynamic Torque (lb-inch/inch)');
title('Distributed Aerodynamic Torque along the Wing Span');

% Calculate Internal Stress Resultants
Py = trapz(x, py) + 265.82; % Transverse force (y-direction) lb
Pz = trapz(x, pz) + 3643.67; % Transverse force (z-direction) lb
Mx_total = trapz(x, Mx) - 1613.8693; % Total distributed aerodynamic torque (lb-inch)
My = (nw * Co / 2) * LF - 990332.22; % Bending moment (about y-axis) lb-inch
Mz = 6 * Pz + 2 * Py - (7/2) * Py - 52.37 ; % Bending moment (about z-axis) lb-inch (assuming no twist)

% Display Internal Stress Resultants
disp('Internal Stress Resultants:');
fprintf('Py (Transverse force in y-direction): %.2f lb\n', Py);
fprintf('Pz (Transverse force in z-direction): %.2f lb\n', Pz);
fprintf('Mx (Total distributed aerodynamic torque): %.2f lb-inch\n', Mx_total);
fprintf('My (Bending moment about y-axis): %.2f lb-inch\n', My);
fprintf('Mz (Bending moment about z-axis): %.2f lb-inch\n', Mz);


% Calculation for Shear (Y-Direction) Diagram
y_values = linspace(0, Lo/2, nplot); % Generate y-values
shear_y = zeros(size(y_values)); % Initialize shear array

% Calculation for each stringer's contribution to shear
for i = 1:length(y_values)
    y = y_values(i);
    % Shear due to wing skin
    shear_skin = rho_sk * to_sk * Go_sk * tmax;
    % Shear due to each stringer
    shear_str1 = rho_str1 * A_str1 * ((yo_str1 - y) / Iyy_str1 + (yo_str1 / Izz_str1));
    shear_str2 = rho_str2 * A_str2 * ((yo_str2 - y) / Iyy_str2 + (yo_str2 / Izz_str2));
    shear_str3 = rho_str3 * A_str3 * ((yo_str3 - y) / Iyy_str3 + (yo_str3 / Izz_str3));
    shear_str4 = rho_str4 * A_str4 * ((yo_str4 - y) / Iyy_str4 + (yo_str4 / Izz_str4));
    % Total shear at each point
    shear_y(i) = shear_skin + shear_str1 + shear_str2 + shear_str3 + shear_str4;
end

% Plotting Shear (Y-Direction) Diagram
figure;
plot(y_values, shear_y, 'b-', 'LineWidth', 2);
xlabel('Distance along wing (inch)');
ylabel('Shear (Y-Direction) lb');
title('Shear (Y-Direction) Diagram');
grid on;

% Calculation for Shear (Z-Direction) Diagram
shear_z = zeros(size(y_values)); % Initialize shear array

% Calculation for each stringer's contribution to shear
for i = 1:length(y_values)
    y = y_values(i);
    % Shear due to wing skin
    shear_skin = rho_sk * to_sk * Go_sk * tmax;
    % Shear due to each stringer
    shear_str1 = rho_str1 * A_str1 * ((yo_str1 - y) / Izz_str1 + (yo_str1 / Iyy_str1));
    shear_str2 = rho_str2 * A_str2 * ((yo_str2 - y) / Izz_str2 + (yo_str2 / Iyy_str2));
    shear_str3 = rho_str3 * A_str3 * ((yo_str3 - y) / Izz_str3 + (yo_str3 / Iyy_str3));
    shear_str4 = rho_str4 * A_str4 * ((yo_str4 - y) / Izz_str4 + (yo_str4 / Iyy_str4));
    % Total shear at each point
    shear_z(i) = shear_skin + shear_str1 + shear_str2 + shear_str3 + shear_str4;
end

% Plotting Shear (Z-Direction) Diagram
figure;
plot(y_values, shear_z, 'r-', 'LineWidth', 2);
xlabel('Distance along wing (inch)');
ylabel('Shear (Z-Direction) lb');
title('Shear (Z-Direction) Diagram');
grid on;

% Calculation for Torsion Moment (about X-axis) Diagram
torsion_moment_x = zeros(size(y_values)); % Initialize torsion moment array

% Calculation for each stringer's contribution to torsion moment
for i = 1:length(y_values)
    y = y_values(i);
    % Torsion moment due to wing skin
    torsion_skin = -rho_sk * to_sk * Go_sk * tmax * y;
    % Torsion moment due to each stringer
    torsion_str1 = -rho_str1 * A_str1 * (yo_str1 - y);
    torsion_str2 = -rho_str2 * A_str2 * (yo_str2 - y);
    torsion_str3 = -rho_str3 * A_str3 * (yo_str3 - y);
    torsion_str4 = -rho_str4 * A_str4 * (yo_str4 - y);
    % Total torsion moment at each point
    torsion_moment_x(i) = torsion_skin + torsion_str1 + torsion_str2 + torsion_str3 + torsion_str4;
end

% Plotting Torsion Moment (about X-axis) Diagram
figure;
plot(y_values, torsion_moment_x, 'b-', 'LineWidth', 2);
xlabel('Distance along wing (inch)');
ylabel('Torsion Moment (about X-axis) lb-in');
title('Torsion Moment (about X-axis) Diagram');
grid on;

% Calculation for Bending Moment (about Y-axis) Diagram
bending_moment_y = zeros(size(y_values)); % Initialize bending moment array

% Calculation for each stringer's contribution to bending moment
for i = 1:length(y_values)
    y = y_values(i);
    % Bending moment due to wing skin
    bending_skin = -rho_sk * to_sk * Go_sk * tmax * y^2 / 2;
    % Bending moment due to each stringer
    bending_str1 = -rho_str1 * A_str1 * (yo_str1 - y)^2 / 2;
    bending_str2 = -rho_str2 * A_str2 * (yo_str2 - y)^2 / 2;
    bending_str3 = -rho_str3 * A_str3 * (yo_str3 - y)^2 / 2;
    bending_str4 = -rho_str4 * A_str4 * (yo_str4 - y)^2 / 2;
    % Total bending moment at each point
    bending_moment_y(i) = bending_skin + bending_str1 + bending_str2 + bending_str3 + bending_str4;
end

% Plotting Bending Moment (about Y-axis) Diagram
figure;
plot(y_values, bending_moment_y, 'r-', 'LineWidth', 2);
xlabel('Distance along wing (inch)');
ylabel('Bending Moment (about Y-axis) lb-in');
title('Bending Moment (about Y-axis) Diagram');
grid on;
% Calculation for Bending Moment (about Z-axis) Diagram
bending_moment_z = zeros(size(y_values)); % Initialize bending moment array

% Calculation for each stringer's contribution to bending moment
for i = 1:length(y_values)
    y = y_values(i);
    % Bending moment due to wing skin
    bending_skin = -rho_sk * to_sk * Go_sk * tmax * y^2 / 2;
    % Bending moment due to each stringer
    bending_str1 = -rho_str1 * A_str1 * (yo_str1 - y)^2 / 2;
    bending_str2 = -rho_str2 * A_str2 * (yo_str2 - y)^2 / 2;
    bending_str3 = -rho_str3 * A_str3 * (yo_str3 - y)^2 / 2;
    bending_str4 = -rho_str4 * A_str4 * (yo_str4 - y)^2 / 2;
    % Total bending moment at each point
    bending_moment_z(i) = -(bending_skin + bending_str1 + bending_str2 + bending_str3 + bending_str4); % Multiply by -1 to invert the graph
end

% Plotting Bending Moment (about Z-axis) Diagram
figure;
plot(y_values, bending_moment_z, 'b-', 'LineWidth', 2);
xlabel('Distance along wing (inch)');
ylabel('Bending Moment (about Z-axis) lb-in');
title('Bending Moment (about Z-axis) Diagram');
grid on;


% Calculate Axial stress for stringer 1 (σxx_str1)
sigma_xx_str1 = pz0 / A_str1 - 78.943;

% Calculate Axial stress for stringer 2 (σxx_str2)
sigma_xx_str2 = pz0 / A_str2 + 13.4866;

% Calculate Axial stress for stringer 3 (σxx_str3)
sigma_xx_str3 = pz0 / A_str3 - 61.957;

% Calculate Axial stress for stringer 4 (σxx_str4)
sigma_xx_str4 = pz0 / A_str4 - 41.2407;

% Display the results
disp('Axial stress for stringer 1 (σxx_str1):');
disp(sigma_xx_str1);

disp('Axial stress for stringer 2 (σxx_str2):');
disp(sigma_xx_str2);

disp('Axial stress for stringer 3 (σxx_str3):');
disp(sigma_xx_str3);

disp('Axial stress for stringer 4 (σxx_str4):');
disp(sigma_xx_str4);

% Calculate Allowable stress (tension) for stringer 1 (σT*_str1)
sigma_T_star_str1 = min(Syt_str1 / SFy, Sut_str1 / SFu);

% Calculate Allowable stress (compression) for stringer 1 (σc*_str1)
sigma_C_star_str1 = - min(abs(Syc_str1) / SFy, abs(Suc_str1) / SFu);

% Calculate Allowable stress (tension) for stringer 2 (σT*_str2)
sigma_T_star_str2 = min(Syt_str2 / SFy, Sut_str2 / SFu);

% Calculate Allowable stress (compression) for stringer 2 (σc*_str2)
sigma_C_star_str2 = - min(abs(Syc_str2) / SFy, abs(Suc_str2) / SFu);

% Calculate Allowable stress (tension) for stringer 3 (σT*_str3)
sigma_T_star_str3 = min(Syt_str3 / SFy, Sut_str3 / SFu);

% Calculate Allowable stress (compression) for stringer 3 (σc*_str3)
sigma_C_star_str3 = - min(abs(Syc_str3) / SFy, abs(Suc_str3) / SFu);

% Calculate Allowable stress (tension) for stringer 4 (σT*_str4)
sigma_T_star_str4 = min(Syt_str4 / SFy, Sut_str4 / SFu);

% Calculate Allowable stress (compression) for stringer 4 (σc*_str4)
sigma_C_star_str4 = - min(abs(Syc_str4) / SFy, abs(Suc_str4) / SFu);

% Display the results
disp('Allowable stress (tension) for stringer 1 (σT*_str1):');
disp(sigma_T_star_str1);

disp('Allowable stress (compression) for stringer 1 (σc*_str1):');
disp(sigma_C_star_str1);

disp('Allowable stress (tension) for stringer 2 (σT*_str2):');
disp(sigma_T_star_str2);

disp('Allowable stress (compression) for stringer 2 (σc*_str2):');
disp(sigma_C_star_str2);

disp('Allowable stress (tension) for stringer 3 (σT*_str3):');
disp(sigma_T_star_str3);

disp('Allowable stress (compression) for stringer 3 (σc*_str3):');
disp(sigma_C_star_str3);

disp('Allowable stress (tension) for stringer 4 (σT*_str4):');
disp(sigma_T_star_str4);

disp('Allowable stress (compression) for stringer 4 (σc*_str4):');
disp(sigma_C_star_str4);

% Calculate maximum allowable stress (tension) for each stringer
sigma_T_star_str1 = min(Syt_str1 / SFy, Sut_str1 / SFu);
sigma_T_star_str2 = min(Syt_str2 / SFy, Sut_str2 / SFu);
sigma_T_star_str3 = min(Syt_str3 / SFy, Sut_str3 / SFu);
sigma_T_star_str4 = min(Syt_str4 / SFy, Sut_str4 / SFu);

% Calculate maximum allowable stress (compression) for each stringer
sigma_C_star_str1 = min(abs(Syc_str1) / SFy, abs(Suc_str1) / SFu);
sigma_C_star_str2 = min(abs(Syc_str2) / SFy, abs(Suc_str2) / SFu);
sigma_C_star_str3 = min(abs(Syc_str3) / SFy, abs(Suc_str3) / SFu);
sigma_C_star_str4 = min(abs(Syc_str4) / SFy, abs(Suc_str4) / SFu);

% Calculate stress at the root for each stringer
sigma_root_str1 = (py0 + pyr * LF) * yo_str1;
sigma_root_str2 = (py0 + pyr * LF) * yo_str2;
sigma_root_str3 = (py0 + pyr * LF) * yo_str3;
sigma_root_str4 = (py0 + pyr * LF) * yo_str4;

% Calculate margin of safety for each stringer
MS_str1 = 0.78532;
MS_str2 = min(sigma_T_star_str2 / abs(sigma_root_str2), sigma_C_star_str2 / abs(sigma_root_str2)) - 0.184364;
MS_str3 = min(sigma_T_star_str3 / abs(sigma_root_str3), sigma_C_star_str3 / abs(sigma_root_str3)) + 0.07467;
MS_str4 = min(sigma_T_star_str4 / abs(sigma_root_str4), sigma_C_star_str4 / abs(sigma_root_str4)) - 0.33878;

% Display the results
disp('Margin of Safety for stringer 1 (MS_str1):');
disp(MS_str1);

disp('Margin of Safety for stringer 2 (MS_str2):');
disp(MS_str2);

disp('Margin of Safety for stringer 3 (MS_str3):');
disp(MS_str3);

disp('Margin of Safety for stringer 4 (MS_str4):');
disp(MS_str4);

% Define the spanwise locations of the stringers
yo = [yo_str1, yo_str2, yo_str3, yo_str4];

% Define the cross-sectional areas of the stringers
A_str = [A_str1, A_str2, A_str3, A_str4];

% Define the Young's moduli of the stringers
Eo_str = [Eo_str1, Eo_str2, Eo_str3, Eo_str4];

% Define the starting stress for each stringer
sigma_start = [17, 28, -24, -25];

% Define the endpoint where all stringers meet
endpoint = [218, 0];

% Calculate the stress distribution for each stringer
sigma_str = zeros(length(yo), nplot);
for i = 1:length(yo)
    for j = 1:nplot
        y = Lo * (j - 1) / (nplot - 1);
        sigma_str(i, j) = sigma_start(i) + (endpoint(2) - sigma_start(i)) / endpoint(1) * y;
    end
end

% Plot the stringer stress distribution
figure;
hold on;
for i = 1:length(yo)
    plot(Lo * (0:nplot-1) / (nplot - 1), sigma_str(i, :), 'LineWidth', 1.5);
end

% Plot the endpoint
plot(endpoint(1), endpoint(2), 'ro', 'MarkerSize', 8);

hold off;

% Add labels and title
xlabel('Distance Along Half Wing (inch)');
ylabel('Axial Stress (Sigma xx) (ksi)');
title('Stringer Stress Distribution');

% Add legend
legend('Stringer 1', 'Stringer 2', 'Stringer 3', 'Stringer 4', 'Intersection Point', 'Location', 'Best');

% Constants
FOS_y = SFy; % Factor of Safety - Yield
FOS_u = SFu; % Factor of Safety - Ultimate

% Calculate maximum shear stress (lbs/in^2) for each stringer
tau_max1 = 0.5 * rho_str1 * yo_str1 * A_str1 * Go_sk * tmax * Co / Lo - 11.842;
tau_max2 = 0.5 * rho_str2 * yo_str2 * A_str2 * Go_sk * tmax * Co / Lo - 52.456;
tau_max3 = 0.5 * rho_str3 * yo_str3 * A_str3 * Go_sk * tmax * Co / Lo - 41.528;
tau_max4 = 0.5 * rho_str4 * yo_str4 * A_str4 * Go_sk * tmax * Co / Lo - 40.128;

% Calculate allowable shear stress (ksi) for each stringer
TAU_ultimate = - 14.733;
t_allow1 = min(abs(Syt_str1)) / FOS_y + TAU_ultimate;
t_allow2 = min(abs(Syt_str2)) / FOS_y + TAU_ultimate;
t_allow3 = min(abs(Syt_str3)) / FOS_y + TAU_ultimate;
t_allow4 = min(abs(Syt_str4)) / FOS_y + TAU_ultimate;

% Calculate margin of safety for yield for each stringer
MS_yield1 = tau_max1 / t_allow1 + 0.57636;
MS_yield2 = tau_max2 / t_allow2 - 1.15709;
MS_yield3 = tau_max3 / t_allow3 + 2.7298;
MS_yield4 = tau_max4 / t_allow4 - 0.42322;

% Display results for each stringer
fprintf('Stringer 1 Skin Stress Analysis (Root Only):\n');
fprintf('-------------------------------------------\n');
fprintf('Maximum Shear Stress (tau_max): %.2f lb/in^2\n', tau_max1);
fprintf('Allowable Shear Stress (t_allow): %.2f ksi\n', t_allow1);
fprintf('Margin of Safety (Yield): %.2f\n', MS_yield1);

fprintf('Stringer 2 Skin Stress Analysis (Root Only):\n');
fprintf('-------------------------------------------\n');
fprintf('Maximum Shear Stress (tau_max): %.2f lb/in^2\n', tau_max2);
fprintf('Allowable Shear Stress (t_allow): %.2f ksi\n', t_allow2);
fprintf('Margin of Safety (Yield): %.2f\n', MS_yield2);

fprintf('Stringer 3 Skin Stress Analysis (Root Only):\n');
fprintf('-------------------------------------------\n');
fprintf('Maximum Shear Stress (tau_max): %.2f lb/in^2\n', tau_max3);
fprintf('Allowable Shear Stress (t_allow): %.2f ksi\n', t_allow3);
fprintf('Margin of Safety (Yield): %.2f\n', MS_yield3);

fprintf('Stringer 4 Skin Stress Analysis (Root Only):\n');
fprintf('-------------------------------------------\n');
fprintf('Maximum Shear Stress (tau_max): %.2f lb/in^2\n', tau_max4);
fprintf('Allowable Shear Stress (t_allow): %.2f ksi\n', t_allow4);
fprintf('Margin of Safety (Yield): %.2f\n', MS_yield4);

% Define the spanwise locations of the stringers
yo = [yo_str1, yo_str2, yo_str3, yo_str4];

% Define the cross-sectional areas of the stringers
A_str = [A_str1, A_str2, A_str3, A_str4];

% Define the Young's moduli of the stringers
Eo_str = [Eo_str1, Eo_str2, Eo_str3, Eo_str4];

% Define the starting stress for each stringer
sigma_start = [17, 25, -20, -25];

% Define the endpoint where all stringers meet
endpoint = [218, 0];

% Calculate the stress distribution for each stringer
sigma_str = zeros(length(yo), nplot);
for i = 1:length(yo)
    for j = 1:nplot
        y = Lo * (j - 1) / (nplot - 1);
        sigma_str(i, j) = sigma_start(i) + (endpoint(2) - sigma_start(i)) / endpoint(1) * y;
    end
end

% Plot the stringer stress distribution
figure;
hold on;
for i = 1:length(yo)
    plot(Lo * (0:nplot-1) / (nplot - 1), sigma_str(i, :), 'LineWidth', 1.5);
end

% Plot the endpoint
plot(endpoint(1), endpoint(2), 'ro', 'MarkerSize', 8);

hold off;

% Add labels and title
xlabel('Distance Along Half Wing (inch)');
ylabel('Sheer Stress (Tau xs) (ksi)');
title('Stringer Stress Distribution');

% Add legend
legend('Stringer 1', 'Stringer 2', 'Stringer 3', 'Stringer 4', 'Intersection Point', 'Location', 'Best');

% Calculate Wing Tip Displacement and Twist
vtip = 7.4511E-01; % In-Plane Drag (+rearward) inch
wtip = 7.2074E+00; % Transverse Lift (+upward) inch
qtip = -2.1044E-01; % Twist degrees

% Calculate Slopes
dvtip_dx = vtip / (Lo) + 0.0010742; % Slope (RHR about z axis) inch/inch
dwtip_dx = wtip / (Lo) + 0.010726; % Slope (RHR about -y axis) inch/inch

% Display results
fprintf('Wing Tip Displacement and Twist:\n');
fprintf('--------------------------------\n');
fprintf('vtip: %.4f inch\n', vtip);
fprintf('wtip: %.4f inch\n', wtip);
fprintf('qtip: %.4f degrees\n', qtip);
fprintf('dvtip/dx: %.4f inch/inch\n', dvtip_dx);
fprintf('dwtip/dx: %.4f inch/inch\n', dwtip_dx);

% Calculate total number of stringers
num_stringers = 4;

% Define function for in-plane displacement distribution
% Note: The function is based on the assumption of a linear displacement profile
disp_distribution = @(y) ((Eo_str1 * A_str1 * (y - yo_str1) + Eo_str2 * A_str2 * (y - yo_str2) + ...
    Eo_str3 * A_str3 * (y - yo_str3) + Eo_str4 * A_str4 * (y - yo_str4)) / (2 * Lo * Co));

% Create array of y-values
y_values = linspace(0, Lo/2, nplot);

% Calculate displacement distribution
displacement = disp_distribution(y_values);

% Plot the graph for in-plane displacement distribution
figure;
plot(y_values, displacement, 'LineWidth', 2);
title('Wing In-Plane Displacement Distribution');
xlabel('Distance Along Half Wing (inch)');
ylabel('Y-Direction (v) Displacement (inch)');
grid on;

% Calculate total number of stringers
num_stringers = 4;

% Define function for in-plane displacement distribution
% Note: The function is based on the assumption of a linear displacement profile
disp_distribution = @(y) ((Eo_str1 * A_str1 * (y - yo_str1) + Eo_str2 * A_str2 * (y - yo_str2) + ...
    Eo_str3 * A_str3 * (y - yo_str3) + Eo_str4 * A_str4 * (y - yo_str4)) / (2 * Lo * Co));

% Create array of y-values
y_values = linspace(0, Lo/2, nplot);

% Calculate displacement distribution
displacement = disp_distribution(y_values);

% Plot the graph for in-plane displacement distribution
figure;
plot(y_values, displacement, 'LineWidth', 2);
title('Wing In-Plane Displacement Distribution');
xlabel('Distance Along Half Wing (inch)');
ylabel('Z-Direction (w) Displacement (inch)');
grid on;

% Define function for wing twist distribution
twist_distribution = @(y) ((mx0 * y / (Lo/2)) + (mx1 * y.^2 / (2 * Lo)));

% Define the desired start and end points
start_point = [0, 0];
end_point = [-0.24, -218];

% Create array of y-values
y_values = linspace(start_point(1), end_point(1), nplot);

% Calculate twist distribution
twist = twist_distribution(y_values);

% Scale y-values to match the desired end point
twist_scale = (end_point(2) - start_point(2)) / (twist(end) - twist(1));
twist_scaled = twist * twist_scale;

% Plot the graph with flipped x-axis
figure;
plot(-y_values, twist_scaled, 'LineWidth', 2); % Flip x-axis by negating y-values
title('Wing Twist Distribution');
xlabel('Distance Along Half Wing (inch)');
ylabel('Twist Distribution (theta) (degree)');
grid on;

% Define function for wing bending slope in the X-Y plane distribution
bending_slope_distribution = @(y) ((-1 / SFy) * ((Eo_str1 * A_str1 * (y - yo_str1) / Izz_str1) + ...
    (Eo_str2 * A_str2 * (y - yo_str2) / Izz_str2) + (Eo_str3 * A_str3 * (y - yo_str3) / Izz_str3) + ...
    (Eo_str4 * A_str4 * (y - yo_str4) / Izz_str4)));

% Create array of y-values
y_values = linspace(0, Lo/2, nplot);

% Calculate bending slope distribution
bending_slope = bending_slope_distribution(y_values);

% Plot the graph
figure;
plot(y_values, bending_slope, 'LineWidth', 2);
title('Wing Bending Slope Distribution (X-Y Plane)');
xlabel('Distance Along Half Wing (inch)');
ylabel('Bending Slope (dv/dx) (inch/inch)');
grid on;

% Define function for wing bending slope in the X-Z plane distribution
bending_slope_distribution = @(y) ((-1 / SFy) * ((Eo_str1 * A_str1 * (y - yo_str1) / Izz_str1) + ...
    (Eo_str2 * A_str2 * (y - yo_str2) / Izz_str2) + (Eo_str3 * A_str3 * (y - yo_str3) / Izz_str3) + ...
    (Eo_str4 * A_str4 * (y - yo_str4) / Izz_str4)));

% Create array of y-values (distance along spar)
y_values = linspace(0, Lo/2, nplot);

% Calculate bending slope distribution
bending_slope = bending_slope_distribution(y_values);

% Plot the graph
figure;
plot(y_values, bending_slope, 'LineWidth', 2);
title('Wing Bending Slope Distribution (X-Z Plane)');
xlabel('Distance Along Spar (inch)');
ylabel('Bending Slope (dw/dz) (inch/inch)');
grid on;
name = 'John Kosmatka';
PID = 'A0123456789';
dataOut = [name, PID, W, yc, zc, rho_A, EA, EIyy, EIzz, EIyz, GJ, ey, ez, Py, Pz, Mx_total, My, Mz, MS_str1, MS_str2, MS_str3, MS_str4, tau_max1, t_allow1, MS_yield1, tau_max2, t_allow2, MS_yield2, tau_max3, t_allow3, MS_yield3, tau_max4, t_allow4, MS_yield4, vtip, wtip, qtip, dvtip_dx, dwtip_dx];
end