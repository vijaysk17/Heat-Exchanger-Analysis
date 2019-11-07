function value = func_const_temp_subsonic(x,M)
    global gamma
    global T_0_1
    global T_wall
    global f
    global choke_check
    
    M = sqrt(M);
    T_0 = T_wall-((T_wall-T_0_1)/exp(2*f*x));                            %Stagnant temperature
    F_T_0 = (M.^2).*(1+(gamma*(M.^2))).*(1+((gamma-1)*(M.^2)/2))./(1-(M.^2));   %Temperature coefficient of differential equation
    F_f = gamma*(M.^4).*(1+((gamma-1)*(M.^2)/2))./(1-(M.^2));                 %Friction coefficient of differential equation
    value = 2*f*((F_T_0*(T_wall-T_0)/T_0)+(2*F_f));                        %derivative of M^2 w.r.t T_0
    choke_check = (M>=1);    
end