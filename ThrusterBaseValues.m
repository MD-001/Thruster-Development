%% Base Thruster Values

%% Temperature

Ti = Tcns / (1 + ((gamma_i - 1) * (Mi^2)) / 2);

Tt = Tcns * ((Pt / Pcns)^((gamma_t - 1) / gamma_t));

Tx = Tcns * ((Px / Pcns)^((gamma - 1) / gamma));

Te = Tcns * ((Pe / Pcns)^((gamma_e - 1) / gamma_e));

%% Sonic Velocity

ai = (gamma_i * R * Ti)^0.5;

at = (gamma_t * R * Tt)^0.5;

ax = (gamma * R * Tx)^0.5;

ae = (gamma_e * R * Te)^0.5;

%% Velocity

vi = Mi*ai;

vt = Mt*at;

vx = (((2*gamma/(gamma-1))*R*Tcns)*...
    (1-(Px/Pcns)^((gamma - 1) / gamma)))^0.5;

ve = (((2*gamma_e/(gamma_e-1))*R*Tcns)*...
    (1-(Pe/Pcns)^((gamma_e-1)/gamma_e)))^0.5 ;


%% Mach Number

Mx = vx/ax;
Me = ve/ae;
M = [Mi, Mt, Mx, Me];

%% Stagnation Temperature

T0i = Tcns / (1 + ((gamma_i - 1)/ 2) * (Mi^2));

T0t = Tcns / (1 + ((gamma_t - 1)/ 2) * (Mt^2));

T0e = Tcns / (1 + ((gamma_e - 1)/ 2) * (Me^2));

StagTemp = [T0i, T0t, T0e];

%% Density

rho_inj = Mr*Pcinj/(R_star*Tcns);
rho_i   = Mr*Pi/(R_star*Ti);
rho_t   = Mr*Pt/(R_star*Tt);
rho_x   = Mr*Px/(R_star*Tx);
rho_e   = Mr*Pe/(R_star*Te);

%% Area 

Ai = mdot_propellant/(vi*rho_i); 
At = mdot_propellant/(vt*rho_t);
Ax = Ax_At*At;
Ae = Ae_At*At;
Ac_At = (1/Mi)*sqrt(((1+((gamma_i-1)/2)*Mi^2)/...
    (1+((gamma_i-1)/2)))^((gamma_i+1)/(gamma_i-1)));
Ac = Ac_At*At;

Area = 250*[Ai, At, Ax, Ae];

%% Performance Indicators
% Characteristic Velocity
    C_star = sqrt(gamma * R * Tcns) / ...
        (gamma * sqrt((2 / (gamma + 1)) ^ ((gamma + 1) / (gamma - 1))));
%C_star = 1688.592; %RS-27

% Theoretical value of Thrust Coefficient
    C_F = sqrt((((2 * gamma ^ 2) / (gamma - 1)) * ...
        ((2 / (gamma + 1)) ^ ((gamma + 1) / (gamma - 1))) * ...
        (1 - (Pe / Pcns) ^ ((gamma - 1) / gamma)))) + ...
        (Ae_At) * ((Pe - P0) / (Pcns));
% Typical Thrust Coefficient in space
    C_F_Space = C_F + (Ae_At) * (P0 / Pcns);

% In space Isp Correction in relation to thrust chamber
%    Istc = C_star * C_F / g0;

% Corrected Characteristic Velocity
    n_cstar = 0.975; % Correction factor for C_star
    C_star_bar = n_cstar * C_star;

% Corrected Thrust Coefficient at Sea level
    n_F = 0.98; % Correction factor for C_F
    C_F_bar = n_F * C_F; %1.890;

% Actual Coefficient of thrust in Space
    C_Fs_actual = n_F * C_F_Space;

% Actual Specific Impulse at Sea level
    Is_tc_actual = (C_star_bar * C_F_bar) / g0;

% Actual Specific Impulse in Space
    Is_tcs_actual = (C_star_bar * C_Fs_actual) / g0;

% Actual Correction factor for exhaust velocity
%   nv = Is_tc_actual / Istc;

% Weight flow rate
    Wdot_tc = mdot_propellant * g0;
% Thrust at Sea (N)
    F_actual_sea = Is_tc_actual * Wdot_tc; 
% Thrust in Space (N)
    F_actual_space = Is_tcs_actual * Wdot_tc; 
% m^2
    At_actual = F_actual_sea / (C_F_bar * Pcns); 
% m^2
    Ae_actual = 12 * At_actual; 

%% Display Values
diary ThursterDesignParameters

disp('---------------------------------------');
fprintf('Initial Design Variables (IPA/LOX) @ %.0f kN\n', F/1000);
disp('---------------------------------------');
fprintf('γ:         %.4f \n', gamma);
fprintf('F:         %.0f kN\n', F/1000);
fprintf('Tcns:      %.2f K\n', Tcns);
fprintf('Pcns:      %.2f bar\n', Pcns/100000);
fprintf('Mr:        %.2f kg/mol\n', Mr);
fprintf('Ae/At:     %.2f \n', Ae_At);
fprintf('ρ_fuel:    %.2f kg/m^3\n', rho_f);
fprintf('mdot_p:    %.2f kg/s\n', mdot_propellant);
fprintf('mdot_o:    %.2f kg/s\n', mdot_oxidiser);
fprintf('mdot_f:    %.2f kg/s\n', mdot_fuel);

disp('---------------------------------------');
disp('Pressure');
disp('---------------------------------------');

fprintf('Pi:        %.1f bar\n', Pi/100000);
fprintf('Pt:        %.1f bar\n', Pt/100000);
fprintf('Pe:        %.1f bar\n', Pe/100000);

disp('---------------------------------------');
disp('Temperature');
disp('---------------------------------------');

fprintf('Ti:        %.1f K\n', Ti);
fprintf('Tt:        %.1f K\n', Tt);
fprintf('Te:        %.1f K\n', Te);

disp('---------------------------------------');
disp('Sonic Velocity');
disp('---------------------------------------');

fprintf('ai:        %.1f m/s\n', ai);
fprintf('at:        %.1f m/s\n', at);
fprintf('ae:        %.1f m/s\n', ae);

disp('---------------------------------------');
disp('Flow Velocity');
disp('---------------------------------------');

fprintf('vi:        %.1f m/s\n', vi);
fprintf('vt:        %.1f m/s\n', vt);
fprintf('ve:        %.1f m/s\n', ve);

disp('---------------------------------------');
disp('Density');
disp('---------------------------------------');

fprintf('ρi:        %.3f kg/m^3\n', rho_i);
fprintf('ρt:        %.3f kg/m^3\n', rho_t);
fprintf('ρe:        %.3f kg/m^3\n', rho_e);

disp('---------------------------------------');
disp('Area');
disp('---------------------------------------');

fprintf('Ai:        %.2f mm²\n', Ai*1000000);
fprintf('At:        %.2f mm²\n', At*1000000);
fprintf('Ae:        %.2f mm²\n', Ae*1000000);

disp('---------------------------------------');
disp('Mach Number');
disp('---------------------------------------');

fprintf('Mi:        %.1f \n', Mi);
fprintf('Mt:        %.1f \n', Mt);
fprintf('Me:        %.1f \n', Me);

disp('---------------------------------------');
disp('Performance Indicators');
disp('---------------------------------------');

fprintf('C_star:    %.1f m/s\n', C_star);
fprintf('C_F:       %.3f \n', C_F);

diary off