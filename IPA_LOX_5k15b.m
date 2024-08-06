%----IPA/LOX 5kN 15bar ----------------------------------------------------
%     
%Design Values
gamma = 1.1565;             % FROM NASA CEA         
gamma_i = 1.1498;           % FROM NASA CEA
gamma_t = 1.1565;           % FROM NASA CEA
gamma_e = 1.2062;           % FROM NASA CEA
gamma_array = [gamma_i, gamma_t, gamma_e];
%
F = 5000;                   % (N) Design Thrust 
I_sp = 241.2;               % IPA/LOX @15bar OFR = 1.5 SeaLevel
Tcns = 3125;                % (K) Combustion Temp  FROM NASA CEA 
Pcns = 1500000;             % (Pa) Combustion Pressure FROM NASA CEA 
Mr = 22.953;                % (kg/mol) Molar Mass of Propellant FROM NASA CEA  
OFR = 1.5;                  % FROM NASA CEA 
Ae_At =  3.02;              % Area Ratio FROM NASA CEA
L_star = 0.3;               % (m) Characteristic Length  (IPA/N2O)
%L_star = 1.8;              % (m) Characteristic Length  (EtOH/N2O) 
rho_f = 786;                % (kg/m^3) Density of IPA (Interpolated)         
rho_o = 1141;               % LO2 @90K



%Molar Gas Constant Nm/kg K 
R  = R_star/Mr;             
cp_propellant = (gamma/(gamma - 1))*R; 

%Mass Flow Rate 
%%(kg/s) Airborne Max IPA Flow (0.5kg/s IPA Max)
%%(kg/s) Airborne Max LOX Flow (7.5kg/s LO2 Max)

    %--Thremodynamic Transport Properties of IPA @ 333K
%cp_coolant =  3117;         % (J/kg K) Specific Heat of Coolant    
%rho_coolant = 748;          % (kg/m^3) Density of Coolant 
%mu_coolant = 0.000824;      % (Ns/m^2) Dynamic Viscosity 
%k_coolant = 0.127;          % (W/mK) Thermal Conductivity
%mdot_coolant =  10;         % (kg/s) Coolant Mass Flow Rate 
%Tco =  333;                 % (K) Water Coolant Temperature

    %--Thremodynamic Transport Properties of Water @ 333K
%cp_coolant =  4183;         % (J/kg K) Specific Heat of Coolant     
%rho_coolant = 983;          % (kg/m^3) Density of Coolant  
%mu_coolant = 0.0004666;     % (Ns/m^2) Dynamic Viscosity
%k_coolant = 0.641;          % (W/mK) Thermal Conductivity
%mdot_coolant =  10;         % (kg/s) Coolant Mass Flow Rate           
%Tco =  333;                 % (K) Water Coolant Temperature

    %--Thremodynamic Transport Properties of Water @ 298K
cp_coolant =  4183;         % (J/kg K) Specific Heat of Coolant 
rho_coolant = 997.1;        % (kg/m^3) Density of Coolant  
mu_coolant = 0.00089050;    % (Ns/m^2) Dynamic Viscosity
k_coolant = 0.5948;         % (W/mK) Thermal Conductivity
mdot_coolant =  1;          % (kg/s) Coolant Mass Flow Rate         
Tco =  298;                 % (K) Water Coolant Temperature