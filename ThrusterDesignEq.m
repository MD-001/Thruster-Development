clear
clc
close all
%Thruster Design Eq Script

%% Control Variables
g0 = 9.80665;               %Standard Acceleration due to Gravity (m/s^2)
P0 = 101325;                %Standard Atmospheric Pressure (kPa)
R_star = 8314.460;          %Gas Constant (Nm kmol^-1K^-1)
Mt = 1;                     %Mach Number at Throat 
Mi = 0.4;                   %Mach Number of Gas Flow at Inlet (0.15 - 0.45)
Ax_At = 2;                  %Expansion Ratio at position X-
 
%% Engine Parameters 
%run('RP1_LOX.m');
%run('RS_27.m');
%run('RS_10B-2.m');
%run('LH2_LOX.m');
run('IPA_LOX_5k15b.m');
%run('IPA_LOX_5k30b.m');

%% Material Parameters 
run('INC718.m');
%run('Ti64.m');
%run('AlSi10Mg.m');
%run('GrCop-42.m');

%%          

mdot_propellant = F/(I_sp*g0);
mdot_fuel = mdot_propellant/(1+OFR);              
mdot_oxidiser = mdot_propellant*OFR/(1+OFR);       

%--- Pressure ratio between Injection plane and Combustion Nozzle Inlet ---
Pcinj_Pcns = (1 + gamma_i*(Mi)^2) /...
    (1 + (gamma_i-1)*(Mi^2)/2)^(gamma_i/(gamma_i-1));

%Combustion Injection Pressure (kPa)
Pcinj = Pcinj_Pcns*Pcns;  

%Static Inlet Pressure (kPa)
Pi = Pcns; %Pcinj/(1+gamma*(Mi^2));        

%Static Throat Pressure (kPa)
Pt = Pcns*(2/(gamma_t+1))^(gamma_t/(gamma_t-1));

Pt_Pcns = Pt/Pcns;
%Pe_Pcns = P0/Pcns;                 %Theroretical

%----Finding Pressure at Px between the throat and exit plane-------------

run('MullerMethod.m');
Pressure = [Pcinj, Pi, Pt, Px, Pe];


%% ----T, a, v, M, Rho, A ---------------------------------------------------

run('ThrusterBaseValues.m');

%% ---- Thrust Chamber Design -----------------------------------------------

%At = F/(C_F_bar*Pcns)
%[Lcyl, Lconv, Ln, Ltotal, Rexit, Rn, Rt, Rc, Dt, Dc] = ...
%    plotNozzleGeometry(At, F, C_F_bar, Pcns, Ae_At, Ac_At, rho_inj, L_star);

run('NozzleGeometry.m');

%Radius Upstream of Throat (m)
R_us = 1.5*Rt;     
%Radius Downstream of Throat (m)
R_ds = 0.382*Rt;  
R_mean = 0.5*(R_us+R_ds);

%% ---- Regenerative Cooling -------------------------------------------


%% Step 1: Calculate Convective Heat Transfer Coefficient (hg)

%Combustion Chamber Temperature Design Value (K)
Tcns_bar = Tcns*n_cstar^2;

%Dynamic Viscocity @ Chamber, Throat, Exit
mu_propellant = (1.184e-7)*(Mr^0.5)*(Tcns_bar^0.6);   
mu_comb = (1.184e-7)*(Mr^0.5)*(T0i^0.6);
mu_throat = (1.184e-7)*(Mr^0.5)*(T0t^0.6);
mu_exit = (1.184e-7)*(Mr^0.5)*(T0e^0.6);

%Prandalt Number
Pr = (4*gamma)/(9*gamma-5); 
Pr_comb = (4*gamma_i)/(9*gamma_i-5);
Pr_throat = (4*gamma_t)/(9*gamma_t-5);
Pr_exit = (4*gamma_e)/(9*gamma_e-5);

A_ratios = [Ac_At, At/At, Ae_At];
Pr_array = [Pr_comb, Pr_throat, Pr_exit];
mu_array = [mu_comb, mu_throat, mu_exit];


Twg_Tcns = 0.8;                 %Wall GAS/STAG  Carbon deposits 
%Twg_Tcns = Twg/Tcns_bar;       %No Carbon Deposit

%Temp Dependence on Viscocity Exponent on hot gas side
w = 0.6;                        %diatomic gasses

% Dimensionless Correction Factor
sigma = zeros(length(A_ratios), length(M));

% Bartz Equation
hg = zeros(length(A_ratios), length(M));

% Heat transfer coefficient on hot gas side
hgc = zeros(length(A_ratios), length(M)); % Added line for hgc

% Thermal Resistance Values
Rd = 3.397e-7 * (exp(8.079 - 1.053 ./ A_ratios)); % Calculate Rd based on A_ratios

labels = {'I', 'T', 'E'}; % Labels for hg values

diary ThursterDesignParameters
disp('---------------------------------------');
disp('Heat Transfer Parameters ');
disp('---------------------------------------');
diary off 

% Find each sigma, hg, hgc and Rd values at different area ratios at Ai, At, 
% Ax and Ae#
for j = 1:length(A_ratios)

    gamma = gamma_array(j); % Get the corresponding gamma value
    
    sigma(j) = (((0.5*(Twg_Tcns)*...
        (1+((gamma-1)/2)*M(j)^2) + 0.5)^(0.8-0.2*w))*...
        ((1+((gamma-1)/2)*M(j)^2)^(0.2*w)))^-1;
%double check
    hg(j) = ((0.026 / (Dt^0.2))* ...
        ((mu_array.^0.2 * cp_propellant) / (Pr_array.^0.6)* ...
        (Pcns/(C_star_bar))^0.8* ...
        (Dt/R_mean)^0.1))*((1/A_ratios(j))^0.9)*sigma(j);
 
    if A_ratios(j) < 1
        Rd(j) = 3.397e-7 * (exp(8.079 - 1.053 / A_ratios(j)));
    else
        Rd(j) = 3.397e-7 * (exp(7.5 - 0.4749 / A_ratios(j)));
    end
    
    hgc(j) = 1 / ((1/hg(j)) + Rd(j)); % Calculate hgc value
% Display Calculated Values
    diary ThursterDesignParameters
   
    fprintf('σ(%s):     %.2f\n', labels{j}, sigma(j));
    fprintf('hg(%s):    %.2f W/m²K\n', labels{j}, hg(j));
    fprintf('hgc(%s):   %.2f W/m²K\n', labels{j}, hgc(j));   
    fprintf('Rd(%s):    %.2e m²K/W\n', labels{j}, Rd(j));
    fprintf('-------------------------\n');
    diary off 
end
%%%%

%% Step 2: Calculate Heat flux on Gas-side (q_g)


%Adiabatic Wall Temp (K) 
Taw = Tcns_bar*0.923; %0.923;%r_turbulent;
Taw_i = T0i*0.923;
Taw_t = T0t*0.923;
Taw_e = T0e*0.923;

%r = Pr^0.33
%Taw = Tcns*(1+r*((gamma-1)/2)*Mx^2)/1+((gamma-1)/2)*Mx^2);

%Heat per unit surface area trasferred across the boundary layer (W/m^2)  
q = hgc(2)*(Taw-Twg);
qi = hgc(1)*(Taw_i-Twg);           
qt = hgc(2)*(Taw_t-Twg);
qe = hgc(3)*(Taw_e-Twg);


%% Step 3: Coolant side heat treansfer coefficient (hc)


%Coolant Wall Side Temp (K)
Twc = Twg - q*t_w/k_material ;        
Twc_i = Twg - qi*t_w/k_material; 
Twc_t = Twg - qt*t_w/k_material; 
Twc_e = Twg - qe*t_w/k_material; 

%Heat transfer coefficient on coolant side (W/m²K)
hc = q/(Twc - Tco);         
hc_i = qi/(Twc_i - Tco);
hc_t = qt/(Twc_t - Tco);  
hc_e = qe/(Twc_e - Tco);  
%hc = [hc_i, hc_t, hc_e]

HP = ((qt+qe)/2)*At;
Pco = 4000000;
%sigma_t = (Pco - Pcns/(d/2))t_w + E*lambda*q*t_w/2*(1-v)*k_m + 6*


% The preceding code aims to determine the diamater of the coolant channels
% to used via numerical approximation using either McCarthy-Wolf for H2
% based fuels or Gneliski for general heat transfer.

%Prandlt Number of Coolant
Pr_coolant = mu_coolant*cp_coolant/k_coolant;

%run('McCarthy-Wolf-Correlation.m');
run('Gneliski_Correlation.m');

% Define the initial guess for d
d_guess = 0.001;

% Set the tolerance for convergence
tolerance = 1e-6;

% Set the maximum number of iterations
max_iterations = 10000;

% Perform the Newton-Raphson iteration
d_solution = d_guess;
iteration = 0;
error = inf;

while error > tolerance && iteration < max_iterations
    % Compute the function value and derivative
    f = equation_1(d_solution);
    df = (equation_1(d_solution + tolerance) - f) / tolerance;
 
    % Update d using the Newton-Raphson formula
    d_solution = d_solution - f / df;
 
    % Calculate the error
    error = abs(f);
 
    % Increment the iteration counter
    iteration = iteration + 1;
end

%Check if convergence was achieved
if error <= tolerance
%  
    d = real(d_solution);
%   
else
    fprintf('Failed to converge\n');
end

%Number of Tubes
N = pi()*(Dt + 0.8*(d+ 2*t_w))/(d+2*t_w);
%Coolant Velocity
v_coolant = (mdot_coolant/rho_coolant)*1/(N*0.125*pi()*d^2);
%Reylonds Number
Reynolds = rho_coolant*v_coolant*d/mu_coolant;
%Darcy Friction Factor
fd = (1.82*log(Reynolds)-1.64)^-2;
    
if (Pr_coolant>=0.1) && (Pr_coolant<=1)
    %Nuselt Number
    Nu = 0.02155*Reynolds^0.8018*Pr^0.7095;
elseif (Pr_coolant>1) && (Pr_coolant<=3)
    Nu = 0.1253*Reynolds^0.8413*Pr^0.6179;
elseif (Pr_coolant>3) && (Pr_coolant<=1000)
    Nu = 0.00881*Reynolds^0.8991*Pr_coolant^0.3911;
elseif (Reynolds>=3000) && (Reynolds<=5e6)   
    Nu = 0.125*fd*(Reynolds - 1000)*Pr_coolant/...
        (1+12.7*sqrt(0.125*fd)*(Pr_coolant^(2/3)-1));
else
    disp('Pr Value out of Range') 
end    

%Q = mdot_coolant*cp_coolant*(298 - 555);

%r_laminar = Pr^0.5;
%r_turbulent = Pr^0.33;
%r = [r_laminar, r_turbulent];   %Local Recovery Factor
                   

%% Display Calculated Values
diary ThursterDesignParameters

disp('---------------------------------------');
%fprintf('Taw  :    %.2f K\n', Taw);
fprintf('Taw_i:    %.2f K\n', Taw_i);
fprintf('Taw_t:    %.2f K\n', Taw_t);
fprintf('Taw_e:    %.2f K\n', Taw_e);

disp('---------------------------------------');
%fprintf('q  :   %.2f \n', q);
fprintf('qi:       %.0f W/m²\n', qi);
fprintf('qt:       %.0f W/m²\n', qt);
fprintf('qe:       %.0f W/m²\n', qe);

disp('---------------------------------------');
%fprintf('Twc  :    %.2f K\n', Twc);
fprintf('Twc_i:    %.2f K\n', Twc_i);
fprintf('Twc_t:    %.2f K\n', Twc_t);
fprintf('Twc_e:    %.2f K\n', Twc_e);

disp('---------------------------------------');
%fprintf('hc  :    %.2f W/m²K\n', hc);
fprintf('hc_i:     %.2f W/m²K\n', hc_i);
fprintf('hc_t:     %.2f W/m²K\n', hc_t);
fprintf('hc_e:     %.2f W/m²K\n', hc_e);

disp('---------------------------------------');
disp('Cooling Material Parameters (Inonel 718)');                                                 
disp('---------------------------------------');
disp(['E:        ', sprintf('%0.2g', E), ' kN/mm²']);
disp(['λ:        ', sprintf('%2.1g', lambda), '1/mmK']);
disp(['v:        ', sprintf('%0.2g', v)]);
disp(['t_w:      ', sprintf('%0.4g', t_w*1000), ' mm']);
disp(['k:        ', sprintf('%0.3g', k_material), ' W/mK']);
fprintf('Twg:      %.2f K\n', Twg);

disp('---------------------------------------');
disp('Thermodynamic Transport Properties (H2O)');
disp('---------------------------------------');
disp(['cp_co:    ', sprintf('%0.4g', cp_coolant), ' J/kg K']);
disp(['ρ_co:     ', sprintf('%0.4g', rho_coolant), ' kg/m^3']);
disp(['µ_co:     ', sprintf('%0.3g', mu_coolant), ' Ns/m²']);
disp(['k_co:     ', sprintf('%0.2g', k_coolant), ' W/mK']);
disp(['mdot_co:  ', sprintf('%0.4g', mdot_coolant), ' kg/s']);
disp(['Tco:      ', sprintf('%0.4g', Tco), ' K']);

disp('---------------------------------------');   
disp('Regenerative Cooling Parameters');
disp('---------------------------------------');
fprintf('Diamater of Tubes: %.2f mm\n', d*1000);
fprintf('Number of Tubes:   %.0f \n', N);
fprintf('Coolant Velocity:  %.2f m/s\n', v_coolant);
fprintf('Prandlt Number:    %.2f \n', Pr_coolant);
fprintf('Reynolds:          %.2f \n',Reynolds);
fprintf('Darcy Friction:    %.4f \n',fd);
fprintf('Nuselt Number:     %.2f \n', Nu);


diary off    
    
%% INJECTOR Design

%run('SimplexSwirlInjector.m');
%run('CoaxialSwirl.m');
%run('PintleInjector.m');
run('DoubletImpingingInjector.m')
