%Gneliski Correlation (Generalised Heat Transfer)
%
%--- Base Equations before subsitution-------------------------------------
%
%v_coolant = (mdot_coolant/rho_coolant)/(0.5*N*0.25*pi()*d^2);
%N = pi()*(Dt + 0.8*(d + 2*t_w))/(d+2*t_w);
%Re = rho_coolant*v_coolant*d/mu_coolant; 
%fd = ((1.82*log(Re)-1.64)^-2);
%fd = 0.3164*Re^-0.25)
%Nu = hc*d/k_coolant
%Nu = (0.125*fd*(Re - 1000)*Pr_coolant)/(1+(12.7*sqrt(0.125*fd)*((Pr_coolant^(2/3))-1)));

equation_1 = @(d) hc_t*d/k_coolant - (0.125*((0.3164*(rho_coolant*...
    (((mdot_coolant/rho_coolant)/(0.5*(pi()*(Dt + 0.8*(d + 2*t_w))/...
    (d+2*t_w))*0.25*pi()*d^2)))*d/mu_coolant)^-0.25))*((rho_coolant*...
    (((mdot_coolant/rho_coolant)/(0.5*(pi()*(Dt + 0.8*(d + 2*t_w))/...
    (d+2*t_w))*0.25*pi()*d^2)))*d/mu_coolant) - 1000)*Pr_coolant)/...
    (1+(12.7*sqrt(0.125*((0.3164*(rho_coolant*...
    (((mdot_coolant/rho_coolant)/...
    (0.5*(pi()*(Dt + 0.8*(d + 2*t_w))/(d+2*t_w))*0.25*pi()*d^2)))*...
    d/mu_coolant)^-0.25)))*((Pr_coolant^(2/3))-1)));