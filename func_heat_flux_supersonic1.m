function value = func_heat_flux_supersonic1(x,M)
    global gamma
    global T_0_1
    global T_diff
    global f
    global choke_check
    global lim_duct_length
    global k
    
    M = sqrt(M);
    T_0 = T_0_1+(T_diff*2*f*x);  %Stagnant temperature
    T_wall = T_0+T_diff;
    F_T_0 = (M^2)*(1+(gamma*(M^2)))*(1+((gamma-1)*(M^2)/2))/(1-(M^2)); %Temperature coefficient of differential equation
    F_f = gamma*(M^4)*(1+((gamma-1)*(M^2)/2))/(1-(M^2)); %Friction coefficient of differential equation
    value = 2*f*((F_T_0*T_diff/T_0)+(2*F_f));  %derivative of M^2 w.r.t x
    if M <= 1
        choke_check = 1;
        lim_duct_length(k) =x;
        k=k+1;
    end
end