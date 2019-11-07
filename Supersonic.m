%CODE TO FIND THE PARAMETERS VARIATION IN HEAD EXCHANGER IN DIFFERENT WALL
%CONDTIONS USING GIVEN INPUT PARAMETERS

%SETTING GLOBAL VARIABLES
global gamma                           %Ratio of specific heats for diatomic gas (air)
global R                               %Specific Gas constant for air
global C_p                             %specific Heat capacity at constant pressure for air
global T_0_1                           %Stagnation temperature at inlet
global T_diff                          %Temperature difference for constant heat flux condition
global T_wall                          %Temperature of wall for constant temperature condition
global f                               %Friction factor to include friction effects
global choke_check                     %for checking choke condition
global lim_duct_length                 %length in which choking occurs
global iter_size                       %iteration size for whole analysis
global k                               %dummy variable
global i

fprintf("\nHEAT EXCHANGER ANALYSIS\n");
fprintf("\nSelect the wall condition for the analysis\n");
fprintf("\n1.CONSTANT WALL TEMPERATURE\n2.CONSTANT HEAT FLUX\n");
wall_condition = input("\nYour choice : \n");

%INITIALISATION OF VARIABLES AND INPUT PARAMETERS
gamma = 1.4;
R = 287;
choke_check = 0;
iter_size = 0.005;
k=1;
C_p = gamma*R/(gamma-1);
val = 0; %to check convergence of solution
%Assigning the input parameters for constant wall condition

if wall_condition ==1
    fprintf("\nENTER THE INPUT PARAMETERS FOR CONSTANT WALL TEMPERATURE CONDITION\n");
    duct_length = input("\nEnter the duct length of heat exchanger in number of duct diameters : \n");
    M_1 = input("\nEnter the Mach number of air at the inlet : \n");
    T_wall = input("\nEnter the wall temperature : \n");
    T_diff = 0;                                    %Not used in constant wall temperature case
    f = input("\nEnter the friction factor :\n");
    lim_duct_length = 0;                           %intialisation of variable
    P_0_1 = input("\nEnter the inlet stagnation pressure in Pa :\n");
    T_0_1 = input("\nEnter the inlet stagnation temperature in Kelvin\n");
    %Outlet pressure should be ambient pressure
    P_3 = input("\nEnter the ambient pressure in Pa :\n");

    
    %Assigning the input parameters for constant Heat flux condition
elseif wall_condition ==2
    fprintf("\nENTER THE INPUT PARAMET100ERS FOR CONSTANT HEAT FLUX CONDITION\n");
    duct_length = input("\nEnter the duct length of heat exchanger in number of duct diameters : \n");
    M_1 = input("\nEnter the Mach number of air at the inlet : \n");
    T_wall = 0;                                    %initialization of variable
    T_diff = input("\nEnter the temperature difference between wall and the flow :\n");
    f = input("\nEnter the friction factor :\n");
    lim_duct_length = 0;                           %intialisation of variable
    P_0_1 = input("\nEnter the inlet stagnation pressure in Pa :\n");
    T_0_1 = input("\nEnter the inlet stagnation temperature in Kelvin\n");
    %Outlet pressure should be ambient pressure
    P_3 = input("\nEnter the ambient pressure in Pa :\n");
    

else
    fprintf("INVALID INPUT!!");
    return;
end

%Initialisation of some more variables
T_1 = T_0_1*power(1+((gamma-1)*(M_1^2)/2),-1);
v_1 = M_1*sqrt(gamma*R*T_1);                    %Velocity of flow at inlet
P_1 = P_0_1*power(T_1/T_0_1,gamma/(gamma-1));
rho_1 = P_1/(R*T_1);                            %Density of air at inlet  
P_ratio_input = P_3/P_1;  

x = 0:iter_size:duct_length;           %DUCT LENGTH ARRAY IN TERMS OF DIAMETER

%SUPERSONIC CONSTANT TEMPERATURE CONDITION
if (wall_condition==1),
    
    fprintf("\nSupersonic constant wall temperature condition\n")
    
    M_sqr_init = M_1^2;
    %INTIAL SOLUTION
    options = odeset('Events',@events_supersonic);
    [x,M_sqr] = ode45(@(x,M)func_const_temp_supersonic1(x,M),x,M_sqr_init);
    
    check = choke_check;
    %SOLVING FOR CHOCKING
    if(choke_check == 1)
       y = 0.001:iter_size:lim_duct_length(1);
       for i = 1:1:length(y)
            choke_check = 0;
            i
            p = 0:iter_size:(lim_duct_length(1)-y(i));
            M_2 = sqrt(M_sqr(length(p),1));
            T_0_2 = T_wall-((T_wall-T_0_1)./exp(2*f*(lim_duct_length(1)-y(i))));
            stg_T_ratio_1 = T_0_2/T_0_1;
            T_ratio_1 = stg_T_ratio_1*(1+((gamma-1)*(M_1^2)/2))/(1+((gamma-1)*(M_2^2)/2));    
            P_ratio_1 = (M_1/M_2)*sqrt(T_ratio_1);
            M_2_shock = sqrt((1+((gamma-1)*M_2/2))/((gamma*M_2^2)-((gamma-1)/2)));
            P_ratio_2 = (1+(gamma*(M_2^2)))/(1+(gamma*(M_2_shock^2)));
            M_sqr_init_afterShock = M_2_shock^2;
            z = (lim_duct_length(1)-y(i)):iter_size:duct_length;
            opts = odeset('Events',@events_subsonic);
            [z,M_sqr_afterShock] = ode45(@(z,M)func_const_temp_supersonic2(z,M),z,M_sqr_init_afterShock);
            if(choke_check==0),
            M_3 = sqrt(M_sqr_afterShock(length(z),1));
            T_0_3 = T_wall-((T_wall-T_0_1)./exp(2*f*duct_length));
            stg_T_ratio_3 = T_0_3/T_0_2;
            T_ratio_3 = stg_T_ratio_3*(1+((gamma-1)*(M_2_shock^2)/2))/(1+((gamma-1)*(M_3^2)/2));    
            P_ratio_3 = (M_2_shock/M_3)*sqrt(T_ratio_3);
            P_ratio_cal = P_ratio_1*P_ratio_2*P_ratio_3;
            else 
            P_ratio_cal = 100000; %%Some random large number    
            end 
            val = abs(P_ratio_cal-P_ratio_input);
       if abs(P_ratio_cal-P_ratio_input)<=0.01
            for j= 1:1:length(p)
                M_sqr_1(j) = M_sqr(j);
            end
            for s= 1:1:(length(x)-length(p))
                M_sqr_1(length(p) + s) = M_sqr_afterShock(s);
            end
            M_sqr =(M_sqr_1)';               %Mach no. profile
            break;
       end
       end
   
    end
    M = sqrt(M_sqr);               %Mach no. profile
    Mach_inlet = M(1);
    T_0 = T_wall-((T_wall-T_0_1)./exp(2*f*x));  %stagnation temperature profile
    T_wall_profile = T_wall*ones(length(x),1);         %wall temperature profile
 end

%SUPERSONIC CONSTANT HEAT FLUX CONDITION
if (wall_condition==2)
    
    fprintf("\nSupersonic constant heat flux condition\n")
    M_sqr_init = M_1^2;
    %INTIAL SOLUTION 
    options = odeset('Events',@events_supersonic);
    [x,M_sqr] = ode45(@(x,M)func_heat_flux_supersonic1(x,M),x,M_sqr_init);
    
    check = choke_check;
    %SOLVING FOR CHOCKING
    if(choke_check == 1)
        y = 0.001:iter_size:lim_duct_length(1);
    for i = 1:1:length(y)
        choke_check = 0;
        i
        p = 0:iter_size:(lim_duct_length(1)-y(i));
        M_2 = sqrt(M_sqr(length(p),1));
        T_0_2 = T_0_1+(T_diff*2*f*(lim_duct_length(1)-y(i)));
        stg_T_ratio_1 = T_0_2/T_0_1;
        T_ratio_1 = stg_T_ratio_1*(1+((gamma-1)*(M_1^2)/2))/(1+((gamma-1)*(M_2^2)/2));    
        P_ratio_1 = (M_1/M_2)*sqrt(T_ratio_1);
        M_2_shock = sqrt((1+((gamma-1)*M_2/2))/((gamma*M_2^2)-((gamma-1)/2)));
        P_ratio_2 = (1+(gamma*(M_2^2)))/(1+(gamma*(M_2_shock^2)));
        M_sqr_init_afterShock = M_2_shock^2;
        z = (lim_duct_length(1)-y(i)):iter_size:duct_length;
        opts = odeset('Events',@events_subsonic);
        [z,M_sqr_afterShock] = ode45(@(z,M)func_heat_flux_supersonic2(z,M),z,M_sqr_init_afterShock,opts);
        if(choke_check==0),
        M_3 = sqrt(M_sqr_afterShock(length(z),1));
        T_0_3 = T_0_1+(T_diff*2*f*(duct_length));
        stg_T_ratio_3 = T_0_3/T_0_2;
        T_ratio_3 = stg_T_ratio_3*(1+((gamma-1)*(M_2_shock^2)/2))/(1+((gamma-1)*(M_3^2)/2));    
        P_ratio_3 = (M_2_shock/M_3)*sqrt(T_ratio_3);
        P_ratio_cal = P_ratio_1*P_ratio_2*P_ratio_3;
        else
            P_ratio_cal = 100000;
        end 
        val = abs(P_ratio_cal-P_ratio_input);
      if abs(P_ratio_cal-P_ratio_input)<=0.01
        for j= 1:1:length(p)
             M_sqr_1(j) = M_sqr(j);
        end
        for s= 1:1:(length(x)-length(p))
             M_sqr_1(length(p) + s) = M_sqr_afterShock(s);
        end
        M_sqr = (M_sqr_1)';               %Mach no. profile
        break;
   end
    end
    end
    M = (sqrt(M_sqr)); %Mach no. profile
    Mach_inlet = M(1);
    T_0 = T_0_1+(T_diff*2*f*x);  %stagnation temperature profile
   
end

if (check==1 & val>0.01) %%Checking final extreme case where normal shock shifts to inlet
    M_3 = 0.9999;
    x = duct_length : -iter_size : 0;
    M_sqr_init_choke = M_3^2;
    if(wall_condition==1)
    [x,M_sqr] = ode45(@(x,M)func_const_temp_subsonic(x,M),x,M_sqr_init_choke);    
    end    
    if(wall_condition==2)
    [x,M_sqr] = ode45(@(x,M)func_const_heat_flux_subsonic(x,M),x,M_sqr_init_choke);
    end
    M  = sqrt(M_sqr);
    Mach_inlet = M(length(M));
end

T_inlet = T_0_1*power(1+((gamma-1)*(Mach_inlet^2)/2),-1);
P_inlet = P_0_1*power(T_inlet/T_0_1,gamma/(gamma-1));
rho_inlet = P_inlet/(R*T_inlet);
c_inlet = sqrt(gamma*R*T_inlet);
v_inlet = Mach_inlet*c_inlet;

l = x;

if(wall_condition==1),
    stag_temp = T_wall-((T_wall-T_0_1))./exp(2*f*l);
end    
if(wall_condition==2),
    stag_temp = T_0_1+(T_diff*2*f*l);
end

%%DEFINING THE RATIOS AND ABSOLUTE VALUES OF VARIOUS FLOW PROPERTIES%%
stag_temp_ratio = stag_temp/T_0_1;  
temp_ratio = stag_temp_ratio.*(1+((gamma-1)*(Mach_inlet^2)/2))./(1+((gamma-1)*(M.^2)/2));
pressure_ratio = (Mach_inlet./M).*sqrt(temp_ratio);
velocity_ratio = (M/Mach_inlet).*sqrt(temp_ratio);
density_ratio = pressure_ratio./temp_ratio;
stag_pressure_ratio = pressure_ratio.*power((1+((gamma-1)*(M.^2)/2))./(1+((gamma-1)*(Mach_inlet^2)/2)),gamma/(gamma-1));
entropy = (C_p*log(stag_temp_ratio))-(R*log(stag_pressure_ratio));

temp = T_inlet*temp_ratio;
pressure = P_inlet*pressure_ratio;
density = rho_inlet*density_ratio;
velocity = v_inlet*velocity_ratio;
stag_pressure = P_0_1*stag_pressure_ratio;
c = velocity./M;
enthalpy = C_p*temp;

%%PLOTTING THE VARIOUS QUANTITIES OF INTEREST%%

figure
subplot(2,2,1)
plot(x,M);
xlabel('Distance in duct diameters');
ylabel('Mach Number');
    
subplot(2,2,2)
plot(l,stag_pressure);
xlabel('Distance in duct diameters');
ylabel('Stagnation pressure in Pa');
    
subplot(2,2,3)
plot(l,temp);
xlabel('Distance in duct diameters');
ylabel('Temperature in K');
    
subplot(2,2,4)
plot(l,stag_temp);
xlabel('Distance in duct diameters');
ylabel('Stagnation temperature in Kelvin');

figure
subplot(3,2,1)
plot(l,pressure);
xlabel('Distance in duct diameters');
ylabel('Pressure in Pa');
            
subplot(3,2,2)
plot(l,density);
xlabel('Distance in duct diameters');
ylabel('Density in kg/m^3');
    
subplot(3,2,3)
plot(l,velocity);
xlabel('Distance in duct diameters');
ylabel('Velocity in m/s');  

subplot(3,2,4)
plot(l,c);
xlabel('Distance in duct diameters');
ylabel('Speed of sound in m/s');

subplot(3,2,5)
plot(entropy,enthalpy);
xlabel('Entropy difference (J/kg)');
ylabel('Enthalpy (J/kg)');

subplot(3,2,6)
plot(l,entropy);
xlabel('Distance in duct diameters');
ylabel('Entropy difference (J/kg)');