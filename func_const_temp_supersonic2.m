function value = func_const_temp_supersonic2(x,Mshock)
    global gamma
    global T_0_1
    global T_wall
    global f
    global choke_check
    
    Mshock = sqrt(Mshock);
    T_0 = T_wall-((T_wall-T_0_1)/exp(2*f*x));  %Stagnant temperature
    F_T_0 = (Mshock^2)*(1+(gamma*(Mshock^2)))*(1+((gamma-1)*(Mshock^2)/2))/(1-(Mshock^2)); %Temperature coefficient of differential equation
    F_f = gamma*(Mshock^4)*(1+((gamma-1)*(Mshock^2)/2))/(1-(Mshock^2));  %Friction coefficient of differential equation
    value = 2*f*((F_T_0*(T_wall-T_0)/T_0)+(2*F_f));                %derivative of M^2 w.r.t x
    if(Mshock>=1)
        choke_check = 1;
    end
    
end