
%----- Cooling Channel Material Parameters (Inconel 718) ------------------
%Inconel is designed to not exceed Twg Value, therefore mechanical properties
%are interpolated for <=Twg

%@ 333K Coolant Temp @555K Interpolation
%E = 1.70e11;                %Youngs Modulus (kN/mm²) (Inonel 718)
%lambda = 14.6e-6;           %Thermal Expansion Coeff (1/mmK) 
%v = 0.272;                  %Poisson's Ratio
%t_w  = 0.000329;            %Wall Thickness (Inonel 718 (0.254 - 1.02 mm)
%k_material = 15.5;          %Thermal conductivity (W/mK)
%Twg =  660;                 %Internal Thrust Chamber Wall Design Temp for carbon based fuels (K)  

    %@ 333K Coolant Temp @600K Interpolation (EXCEL DATA)
E = 1.85E+11;               % (kN/mm²) Youngs Modulus  (Inonel 718)
lambda = 1.43E-05;          % (1/mmK) Thermal Expansion Coeff  
v = 0.332;                  % Poisson's Ratio
k_material = 15.7;          % (W/mK) Thermal conductivity 
Twg =  600;                 % (K) Internal Thrust Chamber Wall Design Temp for carbon based fuels 
t_w  = 0.001;               % (m) Wall Thickness (Inonel 718 (0.254 - 1.02 mm)

    %@ 333K Coolant Temp @555K Interpolation (EXCEL DATA)
%E = 1.88E+11;               % (kN/mm²) Youngs Modulus  (Inonel 718)
%lambda = 1.42E-05;          % (1/mmK) Thermal Expansion Coeff  
%v = 0.329;                  % Poisson's Ratio
%t_w  = 0.0010;              % (m) Wall Thickness (Inonel 718 (0.254 - 1.02 mm)
%k_material = 14.9;          % (W/mK) Thermal conductivity 
%Twg =  555;                 % (K) Internal Thrust Chamber Wall Design Temp for carbon based fuels 

    %@ 75K Coolant Temp @833K Interpolation
%E = 204.9;                 %Youngs Modulus (N/m²) (Inonel 718)
%lambda = 13.8e-6;          %Thermal Expansion Coeff (1/mmK) 
%v = 0.274;                 %Poisson's Ratio          
%t_w =  0.0002;             %Wall Thickness (Inonel 718 (0.254 - 1.02 mm)
%k_material = 20;           %Thermal conductivity (W/mK)
%Twg = 833;                 %Internal Thrust Chamber Wall Design Temp for non-carbon based fuels (K)