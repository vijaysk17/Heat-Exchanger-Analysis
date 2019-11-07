function value = func_heat_flux_supersonic2(x,Mshock)
    global gamma
    global T_0_1
    global T_diff
    global f
    
    Mshock = sqrt(Mshock);
    T_0 = T_0_1+(T_diff*2*f*x);  %Stagnant temperature
    T_wall = T_0+T_diff;
    F_T_0 = (Mshock^2)*(1+(gamma*(Mshock^2)))*(1+((gamma-1)*(Mshock^2)/2))/(1-(Mshock^2));  %Temperature coefficient of differential equation
    F_f = gamma*(Mshock^4)*(1+((gamma-1)*(Mshock^2)/2))/(1-(Mshock^2));  %Friction coefficient of differential equation
    value = 2*f*((F_T_0*T_diff/T_0)+(2*F_f));  %derivative of M^2 w.r.t x
    if(Mshock>=1),
        choke_check = 1;
    end
end