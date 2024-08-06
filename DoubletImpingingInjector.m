%Dimensionless Friction Factor (0.5 - 0.92)
Cd = 0.75; 

%Injection Pressure Drop
%v_o = Cd*sqrt(1.379e6/0.5*(rho_o));
%dPi = (1/2)*(rho_o)*(v_o/Cd)^2
%dPi = Pcns*1.2;  %1.379e6;  %        %*Needs to be ammended
%dP_oxidiser = (1/(2*rho_o))*(mdot_oxidiser/Ao*Cd)^2;
%dP_fuel = (1/(2*rho_f))*(mdot_fuel/Af*Cd)^2;
dP_oxidiser = 0.2*Pcns;
dP_fuel = 0.2*Pcns;

%Total Area of Orifice (m)
Ao = mdot_oxidiser/(Cd*sqrt(2*rho_o*dP_oxidiser));
Af = mdot_fuel/(Cd*sqrt(2*rho_f*dP_fuel));

%Number of Orifices (doublet)
%N_orifice = 700;
N_orifice = 20;

%Area of idividual Orifices (m)
ao = (Ao/(2*N_orifice));
af = (Af/(2*N_orifice));
%ao = ((do^2)*pi/4);
%af = ((df^2)*pi/4);

%Design Injector area for Oxidiser and Fuel
%Diameter of Injection Orifice (m)
do = sqrt(4*ao/pi())*1000;
df = sqrt(4*af/pi())*1000;

%Injection Velocity of Oxidiser and Fuel (m/s)
v_oxidiser = mdot_oxidiser/(Ao*rho_o);
v_fuel = mdot_fuel/(Af*rho_f);

%Injection Momentum Ratio
Rm = mdot_oxidiser*v_oxidiser/(mdot_fuel*v_fuel);

a_f = deg2rad(45);
a_o = deg2rad(45);
% 60° and 45°

Beta = (mdot_fuel*v_fuel*sin(a_f)-mdot_oxidiser*v_oxidiser*sin(a_o))/...
        (mdot_fuel*v_fuel*cos(a_f)+mdot_oxidiser*v_oxidiser*cos(a_o));

%   -Beta = Resultant propellants impact at axis of symmetry
%   +Beta = Resultant propellants impact chamber wall

%liquid hypergolic propellants: 0.03491 - 0.08727 rad
%Crogenic propellants: -Beta

disp('---------------------------------------');
disp('Injector Geometry');
disp('---------------------------------------');
disp(['Fuel Pressure Drop:          ', sprintf('%0.3g',dP_fuel/100000), ' bar']);
disp(['Oxidiser Pressure Drop:      ', sprintf('%0.3g',dP_oxidiser/100000), ' bar']);
disp(['Number of Injectors:         ', num2str(N_orifice)]);
disp(['Oxider Injector Diameter:    ', sprintf('%0.3g',do), ' mm']);
disp(['Fuel Injector Diamater:      ', sprintf('%0.3g',df), ' mm' ]);
disp(['Oxidiser Injection Velocity: ', sprintf('%0.4g',v_oxidiser), ' m/s' ]);
disp(['Fuel Injection Velocity:     ', sprintf('%0.4g',v_fuel),' m/s' ]);  
disp(['Injection Momentum ratio:    ', sprintf('%0.4g',Rm), ]); 
disp(['Resultant Impingement Angle: ', sprintf('%0.4g',rad2deg(Beta)),' °' ])    
