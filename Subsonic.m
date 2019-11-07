%CODE TO FIND THE PARAMETERS VARIATION IN HEAT EXCHANGER IN DIFFERENT WALL
%CONDTIONS USING GIVEN INPUT PARAMETERS

%SETTING GLOBAL VARIABLES
global gamma                           %Ratio of specific heats of diatomic gas (air)
global R                               %Specific Gas constant for air
global C_p                             %specific Heat capacity at constant pressure for air
global T_0_1                           %Stagnation temperature at inlet
global T_diff                          %Temperature difference for constant heat flux condition
global T_wall                          %Temperature of wall for constant temperature condition
global f                               %friction factor to include friction effects
global choke_check                     %for checking choke condition
global lim_duct_length                 %length in which choking occurs
global iter_size                       %iteration size for whole analysis
global k                               %dummy variable
global i
global M_iter

fprintf("\nHEAT EXCHANGER ANALYSIS\n");
fprintf("\nSelect the wall condition for the analysis\n");
fprintf("\n1.CONSTANT WALL TEMPERATURE\n2.CONSTANT HEAT FLUX\n");
wall_condition = input("\nYour choice : \n");

%INITIALISATION OF VARIABLES AND INPUT PARAMETERS
gamma = 1.4;
R = 287;
choke_check = 0;
iter_size = 0.001;
M_iter = 0.0001; 
k=1;
C_p = gamma*R/(gamma-1);

%Assigning the input parameters for constant wall condition

M_1 = 0.0001 : M_iter : 0.9999;    
    

if wall_condition ==1
    fprintf("\nENTER THE INPUT PARAMETERS FOR CONSTANT WALL TEMPERATURE CONDITION\n");
    duct_length = input("\nEnter the duct length of heat exchanger in number of duct diameters : \n");
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
    fprintf("\nENTER THE INPUT PARAMETERS FOR CONSTANT HEAT FLUX CONDITION\n");
    duct_length = input("\nEnter the duct length of heat exchanger in number of duct diameters : \n");
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


for i = 1 : length(M_1),

i
T_1(1,i) = T_0_1*power(1+((gamma-1)*(M_1(i)^2)/2),-1); 

P_1(1,i) = P_0_1*power(T_1(1,i)/T_0_1,gamma/(gamma-1));

%Initialisation of some more variables
v_1(1,i) = M_1(i)*sqrt(gamma*R*T_1(1,i));                    %Velocity of flow at inlet
rho_1(1,i) = P_1(1,i)/(R*T_1(1,i));                            %Density of air at inlet                 %Stagnation pressure at inlet
P_ratio_input(1,i) = P_3/P_1(1,i);           %Ratio of pressure at outlet and input

iter_size = duct_length/1000;
x(:,i) = 0:iter_size:duct_length;           %DUCT LENGTH ARRAY IN TERMS OF DIAMETER

M_sqr_init(1,i) = M_1(i)^2;

%SUBSONIC CONSTANT TEMPERATURE CONDITION
if (wall_condition==1)

    %fprintf("\nSubsonic constant wall temperature condition\n")
    T_0(:,i) = T_wall-((T_wall-T_0_1))./exp(2*f*x(:,i));
    stg_T_ratio(:,i) = T_0(:,i)/T_0_1;
    
    %INTIAL SOLUTION
    options = odeset('Events',@events_subsonic);
    [x(:,i),M_sqr(:,i)] = ode45(@(x,M)func_const_temp_subsonic(x,M),x(:,i),M_sqr_init(1,i));
    M(:,i) = sqrt(M_sqr(:,i));
    if(choke_check==0)
    M_1_s(1,i) = M_1(i);
    T_ratio(:,i) = stg_T_ratio(:,i).*(1+((gamma-1)*(M_1_s(1,i)^2)/2))./(1+((gamma-1)*(M(:,i).^2)/2));
    P_ratio(:,i) = (M_1_s(1,i)./M(:,i)).*sqrt(T_ratio(:,i));
    P_ratio_cal(i) = P_ratio(length(x(:,i)),i);
    end
    
    %CHECKING FOR CHOKING
    if(choke_check == 1)
        break;
    end
      %stagnation temperature profile
    %T_wall_profile = T_wall*ones(length(x(:,i)),1);         %wall temperature profile               %Mach no. profile
    
end


%SUBSONIC CONSTANT HEAT FLUX CONDITION
if (wall_condition==2)
        
    T_0(:,i) = T_0_1+(T_diff*2*f*x(:,i));
  
    stg_T_ratio(:,i) = T_0(:,i)/T_0_1;
    
    %INTIAL SOLUTION 
    options = odeset('Events',@events_subsonic);
    [x(:,i),M_sqr(:,i)] = ode45(@(x,M)func_const_heat_flux_subsonic(x,M),x(:,i),M_sqr_init(1,i));
    M(:,i) = sqrt(M_sqr(:,i));
    if(choke_check==0)
    M_1_s(1,i) = M_1(i);
    T_ratio(:,i) = stg_T_ratio(:,i).*(1+((gamma-1)*(M_1_s(1,i)^2)/2))./(1+((gamma-1)*(M(:,i).^2)/2));
    P_ratio(:,i) = (M_1_s(1,i)./M(:,i)).*sqrt(T_ratio(:,i));
    P_ratio_cal(1,i) = P_ratio(length(x(:,i)),i);
    end
    
    %CHECKING FOR CHOKING
    if(choke_check == 1)
     break;
    end
   %stagnation temperature profile
    %T_wall_profile = T_0+T_diff;         %wall temperature profile
    %M(:,i) = sqrt(M_sqr(:,i));               %Mach no. profile
end

convergence(1,i) = abs(P_ratio_cal(1,i)-P_ratio_input(1,i));
end

%%CHECKING IF ANY OF THE SOLUTIONS CONVERGED TO THE RIGHT VALUE##
[val,I] = min(convergence);
Mach = M(:,I);
l = x(:,I);
Mach_inlet = Mach(1,1);

%%CHECKING FOR CHOKED FLOW##
if(I>=length(M_1_s) || val>=0.005) ,
    M_3 = 0.9999;
    x(:,i) = duct_length : -iter_size : 0;
    M_sqr_init_choke = M_3^2;
    if(wall_condition==1)
    [x(:,i),M_sqr(:,i)] = ode45(@(x,M)func_const_temp_subsonic(x,M),x(:,i),M_sqr_init_choke);    
    end    
    if(wall_condition==2)
    [x(:,i),M_sqr(:,i)] = ode45(@(x,M)func_const_heat_flux_subsonic(x,M),x(:,i),M_sqr_init_choke);
    end
    M_1_s(1,i) = sqrt(min(M_sqr(:,i)));
    M(:,i) = sqrt(M_sqr(:,i));
    l = x(:,i);
    Mach = M(:,i);
    Mach_inlet = Mach(length(Mach),1);
end
T_inlet = T_0_1*power(1+((gamma-1)*(Mach_inlet^2)/2),-1);
P_inlet = P_0_1*power(T_inlet/T_0_1,gamma/(gamma-1));
rho_inlet = P_inlet/(R*T_inlet);
c_inlet = sqrt(gamma*R*T_inlet);
v_inlet = Mach_inlet*c_inlet;

if(wall_condition==1),
    stag_temp = T_wall-((T_wall-T_0_1))./exp(2*f*l);
end    
if(wall_condition==2),
    stag_temp = T_0_1+(T_diff*2*f*l);
end

%%DEFINING THE RATIOS AND ABSOLUTE VALUES OF VARIOUS FLOW PROPERTIES%%
stag_temp_ratio = stag_temp/T_0_1;  
temp_ratio = stag_temp_ratio.*(1+((gamma-1)*(Mach_inlet^2)/2))./(1+((gamma-1)*(Mach.^2)/2));
pressure_ratio = (Mach_inlet./Mach).*sqrt(temp_ratio);
velocity_ratio = (Mach/Mach_inlet).*sqrt(temp_ratio);
density_ratio = pressure_ratio./temp_ratio;
stag_pressure_ratio = pressure_ratio.*power((1+((gamma-1)*(Mach.^2)/2))./(1+((gamma-1)*(Mach_inlet^2)/2)),gamma/(gamma-1));
entropy = (C_p*log(stag_temp_ratio))-(R*log(stag_pressure_ratio));

temp = T_inlet*temp_ratio;
pressure = P_inlet*pressure_ratio;
density = rho_inlet*density_ratio;
velocity = v_inlet*velocity_ratio;
stag_pressure = P_0_1*stag_pressure_ratio;
c = velocity./Mach;
enthalpy = C_p*temp;

%%PLOTTING THE VARIOUS QUANTITIES OF INTEREST%%

figure
subplot(2,2,1)
plot(l,Mach);
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
ylabel('Entropy difference (J/kg)');
xlabel('Distance in duct diameters');