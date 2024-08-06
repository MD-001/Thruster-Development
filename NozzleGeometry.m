
% Function to calculate and plot the design parameters for a bell-shaped nozzle
% based on the Rao Nozzle design at 80% fractional length.

%% Constants and Parameters
n = 80;                    % Fractional length (%)
alpha = pi()/12;           % Semi-aperture angle (radians)
theta = pi()/9;            % Convergence angle (radians)
theta_n = 0.4782;          % Nozzle inlet angle (radians)
theta_e = 0.1710;          % Nozzle exit angle (radians)

%% Calculated Variables
At = F / (C_F_bar * Pcns); % Throat area (m^2)

% Throat Dimensions
Rt = sqrt(At / pi());     % Radius of the throat (m)
Dt = 2 * Rt;              % Diameter of the throat (m)

% Upstream Nozzle Radius
Rn = 1.5 * Rt;            % Radius of the nozzle upstream (m)

% Combustion Chamber Dimensions
Rc = sqrt(Ac_At) * Rt;    % Radius of combustion chamber (m)
Dc = sqrt(Ae_At);         % Diameter of combustion chamber (m)
Rexit = sqrt(Ae_At) * Rt; % Radius of nozzle exit (m)

% Volume Calculations
Vc = L_star * At;         % Volume of combustion chamber and convergence (m^3)
Lconv = (Rt * (sqrt(Ac_At) - 1) + Rn * (1 / cos(theta) - 1)) / tan(theta); % Length of convergence (m)
Vconv = (1/3) * pi() * Lconv * (Rc^2 + Rc * Rt + Rt^2); % Volume of convergence (m^3)
Vcyl = Vc - Vconv;        % Volume of cylindrical combustion chamber (m^3)
Lcyl = Vcyl / (Ac_At * At); % Length of cylindrical combustion chamber (m)
Ltotal = Lcyl + Lconv;    % Total length from injector to throat (m)

% Diverging Portion Length
Ln = (n / 100) * (Rt / tan(0.2618)) * (sqrt(Ae_At) - 1 + 1.5 * (1 / cos(0.2618) - 1)); % Length of diverging portion (m)

% Nozzle Points
xN = 0.382 * Rt * sin(alpha); % x-coordinate of nozzle inflection point (m)
yN = Rt * (1 + 0.382 * (1 - cos(alpha))); % y-coordinate of nozzle inflection point (m)
xE = Ln;                     % x-coordinate of nozzle exit (m)
yE = sqrt(Ae_At) * Rt;       % y-coordinate of nozzle exit (m)

xQ = (((yE-tan(theta_e)*xE)-(yN-tan(theta_n)*xN))/...
          (tan(theta_n) - tan(theta_e)));
yQ = (xQ*tan(theta_n) + (yN - xN*tan(theta_n)));

% Solve for nozzle curve polynomial coefficients
A1 = [xN^3, xN^2, xN, 1;
      xE^3, xE^2, xE, 1;
      3 * xN^2, 2 * xN, 1, 0;
      3 * xE^2, 2 * xE, 1, 0];
B1 = [yN; yE; tan(theta_n); tan(theta_e)];
coefficients = A1 \ B1;
a = coefficients(1);
b = coefficients(2);
c = coefficients(3);
d = coefficients(4);

% Convert to millimeters
Rt_mm = Rt * 1000;
xN_mm = xN * 1000;
yN_mm = yN * 1000;
xE_mm = xE * 1000;
yE_mm = yE * 1000;
xQ_mm = xQ * 1000;
yQ_mm = yQ * 1000;
Ln_mm = Ln * 1000;

% Generate nozzle curve points in millimeters
xa = linspace(xN_mm, xE_mm, 100);
ya = a * (xa / 1000).^3 + b * (xa / 1000).^2 + c * (xa / 1000) + d;
ya = ya * 1000;

% Equations for circular nozzle segments
equation1 = @(xb, yb) xb.^2 + (yb - 1.382 * Rt_mm).^2 - (0.382 * Rt_mm)^2; % Downstream circle equation
equation2 = @(xc, yc) xc.^2 + (yc - 2.5 * Rt_mm).^2 - (1.5 * Rt_mm)^2;   % Upstream circle equation

% Generate linearly spaced values for plotting
xC_to_xT = linspace(-Lconv * 1000, 0, 100);
yC_to_yT = linspace(2.5 * Rt_mm - 1.5 * Rt_mm, Rc * 1000, 100);
x_ic = -sqrt((1.5 * Rt_mm)^2 - (Rc * 1000 - 2.5 * Rt_mm)^2); % Initial x position for combustion chamber
xi_to_xC = linspace(x_ic, -Lcyl * 1000, 100);
yi_to_yC = linspace(Rc * 1000, Rc * 1000, 100);

% Plot Nozzle Geometry in mm
[xb, yb] = meshgrid(linspace(0, xN_mm, 100), linspace(Rt_mm, yN_mm, 100));
Zb = equation1(xb, yb);
[xc, yc] = meshgrid(xC_to_xT, yC_to_yT);
Zc = equation2(xc, yc);

figure;
hold on;
plot(xa, ya, 'b', 'LineWidth', 2);           % Plot nozzle curve
contour(xb, yb, Zb, [0 0], 'r', 'LineWidth', 2); % Downstream circle
contour(xc, yc, Zc, [0 0], 'g', 'LineWidth', 2); % Upstream circle
plot(xi_to_xC, yi_to_yC, 'LineWidth', 2);    % Combustion chamber section
hold off;
xlabel('x [mm]');
ylabel('y [mm]');
title('Nozzle Geometry');
grid on;
axis equal; % Set equal scaling for both axes
ylim([0, max([ya, yb(:)', yc(:)', yi_to_yC])]);


%% Display Calculated Values
diary ThursterDesignParameters

disp('---------------------------------------');
disp('Nozzle Points');
disp('---------------------------------------');
disp(['xN:        ', sprintf('%0.5g', xN_mm), ' mm']);
disp(['yN:        ', sprintf('%0.4g', yN_mm), ' mm']);
disp(['xQ:        ', sprintf('%0.4g', xQ_mm), ' mm']);
disp(['yQ:        ', sprintf('%0.5g', yQ_mm), ' mm']);
disp(['xE:        ', sprintf('%0.4g', xE_mm), ' mm']);
disp(['yE:        ', sprintf('%0.4g', yE_mm), ' mm']);

disp('---------------------------------------');
disp('Thruster Geometry');
disp('---------------------------------------');
disp(['Lcyl:      ', sprintf('%0.5g', Lcyl * 1000), ' mm']);
disp(['Lconv:     ', sprintf('%0.4g', Lconv * 1000), ' mm']);
disp(['Ln:        ', sprintf('%0.4g', Ln_mm), ' mm']);
disp(['Ltotal:    ', sprintf('%0.5g', Ltotal * 1000), ' mm']);
disp(['Re:        ', sprintf('%0.4g', Rexit * 1000), ' mm']);
disp(['Rn:        ', sprintf('%0.4g', Rn * 1000), ' mm']);
disp(['Rt:        ', sprintf('%0.4g', Rt_mm), ' mm']);
disp(['Dt:        ', sprintf('%0.4g', Dt * 1000), ' mm']);
disp(['Rc:        ', sprintf('%0.4g', Rc * 1000), ' mm']);

disp('---------------------------------------');
disp('Nozzle Polynomial: ');
disp(['(', num2str(a), ')*x^3 + (', num2str(b), ')*x^2 + (', num2str(c), ')*x + (', num2str(d), ')']);

diary off 
