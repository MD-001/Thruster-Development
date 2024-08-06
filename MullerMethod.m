%This function is used to determine an Unkown presure (Px) at any position along
%the nozzle from a predetermined Ax/At point to subsequently find the exit
%presure (Pe) utlisiling the Muller Method to find the roots through
%quadratic interpolation. Details can be found on page 21 of "Fundemetal
%concepts of Liquid Propellant Rocket Engines"

    % Finding Pressure at Px between the throat and exit plane
    zx = Pt / Pcns;    % Exit plane pressure to 4 s.f.

    % Loop through values of z
    % Calculates the z_value for the isentropic function by using smaller
    % values of z until the f_funtion spits outs the first negative
    z_values = zx:-0.0005:0.0005;
    f_values = (Ax_At^2) - ...
        ((2/(gamma+1))^(2/(gamma-1)) * (1./z_values).^(2/gamma)) ./ ...
        (((gamma+1)/(gamma-1)) * (1 - z_values.^((gamma-1)/gamma)));

    % Plotting z_values versus f_values to find point of sign flip
    %plot(z_values, f_values)
    %xlabel('z_values')
    %ylabel('f_values')
    %title('z_values vs f_values')
    %grid on
    %xlim([0, 0.025])   % Set x-axis limits from 0 to 0.025
    %ylim([-2000, 2000])   % Set y-axis limits from -2000 to 2000
%
    %hold on
    
    % Find the first negative f_value and corresponding z_value
    f2pos = find(f_values < 0, 1);
    f2 = f_values(f2pos);
    z2 = z_values(f2pos);

    % Find the last positive f_value and corresponding z_value
    f1pos = find(f_values > 0, 1, 'last');
    f1 = f_values(f1pos);
    z1 = z_values(f1pos);
    
    %scatter(z2, f2, 'r', 'filled')
    %scatter(z1, f1, 'r', 'filled')
    
    % Calculate the value of z between z2 and z1
    z0 = (z2 + z1) / 2;
    f0 = (Ax_At^2) - ...
        ((2/(gamma+1))^(2/(gamma-1)) * (1/z0)^(2/gamma)) / ...
        (((gamma+1)/(gamma-1)) * (1 - z0^((gamma-1)/gamma)));

    h1 = z1 - z0;
    h2 = z0 - z2;
    k = h2 / h1;

    % Find corresponding quadratic polynomial coefficients via interpolation
    A = (k*f1 - f0*(1+k) + f2) / (k*(h1^2)*(1+k));
    B = (f1 - f0 - (A*(h1^2))) / h1;
    C = f0;

    % Estimate the roots for z and corresponding Px Value
    z = z0 - (2*C) ./ (B + sqrt(B^2 - 4*A*C));
    Px = z * Pcns;

    % Finding Pressure at Pe at the exit plane
    ze = Pt / Pcns;    % Exit plane pressure to 4 s.f.

    % Loop through values of z
    ze_values = ze:-0.0001:0.0001;
    fe_values = (Ae_At^2) - ...
        ((2/(gamma+1))^(2/(gamma-1)) * (1./ze_values).^(2/gamma)) ./ ...
        (((gamma+1)/(gamma-1)) * (1 - ze_values.^((gamma-1)/gamma)));

    % Find the first negative f_value and corresponding z_value
    fe2pos = find(fe_values < 0, 1);
    fe2 = fe_values(fe2pos);
    ze2 = ze_values(fe2pos);

    % Find the last positive f_value and corresponding z_value
    fe1pos = find(fe_values > 0, 1, 'last');
    fe1 = fe_values(fe1pos);
    ze1 = ze_values(fe1pos);

    % Calculate the value of z between z2 and z1
    ze0 = (ze2 + ze1) / 2;
    fe0 = (Ae_At^2) - ...
        ((2/(gamma+1))^(2/(gamma-1)) * (1/ze0)^(2/gamma)) / ...
        (((gamma+1)/(gamma-1)) * (1 - ze0^((gamma-1)/gamma)));

    he1 = ze1 - ze0;
    he2 = ze0 - ze2;
    ke = he2 / he1;

    Ae = (ke*fe1 - fe0*(1+ke) + fe2) / (ke*(he1^2)*(1+ke));
    Be = (fe1 - fe0 - (Ae*(he1^2))) / he1;
    Ce = fe0;

    zei = ze0 - (2*Ce) ./ (Be + sqrt(Be^2 - 4*Ae*Ce));
  
    Pe = zei * Pcns;
     
