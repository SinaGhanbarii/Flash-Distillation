%% Unit Operation 1 project - Flash distillation column
% Mohammad Sina Ghanbari Pakdehi        Student ID. = 98109978    Spring 2022 
%% Part1 - The binary systems data
clear,clc
Water_Coeff =       [5.08354,1663.125,-45.622];            % From 344 K to 373 K
Methanol_Coeff =    [5.20409,1581.341,-33.50];             % From 288.1 K to 356.83 K
Benzene_Coeff =     [4.01814,1203.835,-53.226];            % From 287.7 K to 354.07 K
Toluene_Coeff =     [4.07827,1343.943,-53.773];            % From 308.52 K to 384.66 K
Octane_Coeff =      [4.04867,1355.126,-63.633];            % From 326.08 K to 399.72 K
Heptane_Coeff =     [4.02832,1268.636,-56.199];            % From 299.07 K to 372.43 K
CycloHexane_Coeff = [3.96988,1203.526,-50.287];            % From 293.06 K to 354.73 K
Propane_Coeff =     [4.53678,1149.36,24.906];              % From 277.6 K to 360.8 K
Butane_Coeff =      [4.35576,1175.581,-2.071];             % From 272.66 K to 425 K

Cp_Water_mol = 75.95; Cp_Methanol_mol = 81.11; Cp_Benzene_mol = 133; Cp_Toluene_mol = 103.7;
Cp_Octane_mol = 252.4; Cp_Heptane_mol = 224.8; Cp_CycloHexane_mol = 156.7; Cp_Propane_mol = 73.6; Cp_Butane_mol = 129.7;

MW_Water =18.01528; MW_Methanol = 32.04; MW_Benzene = 78.11; MW_Toluene = 92.14; MW_Octane = 114.23;
MW_Heptane = 100.21; MW_CycloHexane = 84.16; MW_Propane = 44.1; MW_Butane = 58.12;

Cp_Water = (Cp_Water_mol*1000)/MW_Water; Cp_Methanol = (Cp_Methanol_mol*1000)/MW_Methanol; Cp_Benzene = (Cp_Benzene_mol*1000)/MW_Benzene;
Cp_Toluene = (Cp_Toluene_mol*1000)/MW_Toluene; Cp_Octane = (Cp_Octane_mol*1000)/MW_Octane; Cp_Heptane = (Cp_Heptane_mol*1000)/MW_Heptane;
Cp_CycloHexane = (Cp_CycloHexane_mol*1000)/MW_CycloHexane; Cp_Propane = (Cp_Propane_mol*1000)/MW_Propane; Cp_Butane = (Cp_Butane_mol*1000)/MW_Butane;

hfg_Water = 2250.76e3; hfg_Methanol = 1035e3; hfg_Benzene = 30.8e3; hfg_Toluene = 38e3; hfg_Octane = 34.41*MW_Octane*1000; hfg_Heptane = 318e3;
hfg_CycloHexane = 33.1*MW_CycloHexane*1000; hfg_Propane = 428e3; hfg_Butane = 386e3;
%% Part2 - Getting input data from user
disp('Hello and welcome to the flash distillation visulization program. There are several binary systems that showed below:')
number = {1;2;3;4;5;6};
Component1 = {'Water';'Benzene';'Toluene';'Heptane';'Propane';'Toluene'};
Component2 = {'Methanol';'Toluene';'Octane';'CycloHexane';'Butane';'CycloHexane'};
table1 = table(number,Component1,Component2); disp(table1)
command1 = input('Please choose the number of your desired system: ');
command2 = input('What is the unit of feed pressure? bar[1] or atm[2]? ');
if command2 == 1
   feed_p_bar = input('Enter the feed pressure in bar: ');
   feed_p_atm = 0.986923*feed_p_bar;
elseif command2 == 2
    feed_p_atm = input('Enter the feed pressure in atm: ');
    feed_p_bar = 1.01325*feed_p_atm;
else
    disp(command2)
end
command3 = input('What is the unit of the feed rate: kmol/h[1] or kg/h[2]: ');
if command3 == 1
   feed_rate = input('Enter the feed rate in kmol/h: ')*1000; 
elseif command3 == 2
   feed_rate = input('Enter the feed rate in kg/h: '); 
else
   disp(command3)
end
command4 = input('What is the type of feed steam? LPS[1] MPS[2] HPS[3]: ');
if command4 == 1
   hfg = 2108.5e3;                 % Vaporization Enthalpy [=] J/kg
elseif command4 == 2    
    hfg = 1988.1e3;                % Vaporization Enthalpy [=] J/kg                
elseif command4 == 3
    hfg = 1714.1e3;                % Vaporization Enthalpy [=] J/kg
else
    disp(command4)
end
feed_q = input('Enter the feed thermal condition: ');
feed_T = input('Enter the feed temperature in Celcius: ');
feed_T_K = feed_T + 273.15;
feed_Z = input('Enter the feed composition in mol/mol: ');

%% Part3 - Calculating equilibrum x & y & Bubble Temperature 
x_eq = linspace(0,1,1000);          % vector of equilibrum_x
T_Bub = ones(1,length(x_eq));       % preallocated vector of Bubble Temp.
y_eq = ones(1,length(x_eq));        % preallocated vector of equilibrum_y.    
if command1 == 1
    for i = 1:1000
        Bubble_T = @(T) TotalPressure(Water_Coeff(1),Water_Coeff(2),Water_Coeff(3),Methanol_Coeff(1),Methanol_Coeff(2),Methanol_Coeff(3) ... 
            ,T,feed_p_bar,x_eq(i));         
        x0 = 300;                               % Initial guess for Bubble T
        T_Bub(1,i) = fsolve(Bubble_T,x0);       % Calculating Bubble Temp. using itrative calculation
        y_eq(1,i) = x_eq(1,i) * Antoine(Methanol_Coeff(1),Methanol_Coeff(2),Methanol_Coeff(3),T_Bub(1,i))/feed_p_bar;
                                                % Calculating equilibrum y
        clc;
    end
elseif command1 == 2
        for i = 1:1000
        Bubble_T = @(T) TotalPressure(Benzene_Coeff(1),Benzene_Coeff(2),Benzene_Coeff(3),Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3) ... 
            ,T,feed_p_bar,x_eq(i));
        x0 = 300;                               % Initial guess for Bubble T
        T_Bub(i) = fsolve(Bubble_T,x0);         % Calculating Bubble Temp. using itrative calculation
        y_eq(1,i) = x_eq(1,i) * Antoine(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),T_Bub(i)) / feed_p_bar;
                                                % Calculating equilibrum y
        clc;
        end
elseif command1 == 3
        for i = 1:1000
        Bubble_T = @(T) TotalPressure(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),Octane_Coeff(1),Octane_Coeff(2),Octane_Coeff(3) ... 
            ,T,feed_p_bar,x_eq(i));
        x0 = 300;                               % Initial guess for Bubble T
        T_Bub(i) = fsolve(Bubble_T,x0);         % Calculating Bubble Temp. using itrative calculation
        y_eq(1,i) = x_eq(1,i) * Antoine(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),T_Bub(i)) /feed_p_bar;
                                                % Calculating equilibrum y
        clc;
        end
elseif command1 == 4
        for i = 1:1000
        Bubble_T = @(T) TotalPressure(Heptane_Coeff(1),Heptane_Coeff(2),Heptane_Coeff(3),CycloHexane_Coeff(1),CycloHexane_Coeff(2),CycloHexane_Coeff(3) ... 
            ,T,feed_p_bar,x_eq(i));
        x0 = 300;                               % Initial guess for Bubble T
        T_Bub(i) = fsolve(Bubble_T,x0);         % Calculating Bubble Temp. using itrative calculation
        y_eq(1,i) = x_eq(1,i) * Antoine(CycloHexane_Coeff(1),CycloHexane_Coeff(2),CycloHexane_Coeff(3),T_Bub(i))/feed_p_bar;
                                                % Calculating equilibrum y
        clc;
        end
elseif command1 == 5
        for i = 1:1000
        Bubble_T = @(T) TotalPressure(Propane_Coeff(1),Propane_Coeff(2),Propane_Coeff(3),Butane_Coeff(1),Butane_Coeff(2),Butane_Coeff(3) ... 
            ,T,feed_p_bar,x_eq(i));
        x0 = 300;                               % Initial guess for Bubble T
        T_Bub(i) = fsolve(Bubble_T,x0);         % Calculating Bubble Temp. using itrative calculation
        y_eq(1,i) = x_eq(1,i) * Antoine(Propane_Coeff(1),Propane_Coeff(2),Propane_Coeff(3),T_Bub(i)) /feed_p_bar;
                                                % Calculating equilibrum y
        clc;
        end
elseif command1 == 6
        for i = 1:1000
        Bubble_T = @(T) TotalPressure(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),CycloHexane_Coeff(1),CycloHexane_Coeff(2),CycloHexane_Coeff(3) ... 
            ,T,feed_p_bar,x_eq(i));
        x0 = 300;                               % Initial guess for Bubble T
        T_Bub(i) = fsolve(Bubble_T,x0);         % Calculating Bubble Temp. using itrative calculation
        y_eq(1,i) = x_eq(1,i) * Antoine(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),T_Bub(i)) /feed_p_bar;
                                                % Calculating equilibrum y
        clc;
        end
else
    disp(command1)
end

%% Part4 - Calculating the vapor & liquid product amount & composition
% Calculating the value of liquid & vapor product
error = ones(1,length(x_eq));       % preallocated vector of error.
if feed_q >=1 
   L_out = feed_rate;
   V_out = 0;
elseif feed_q <= 0
    L_out = 0;
    V_out = feed_rate;
else
    L_out = feed_q*feed_rate;
    V_out = (1-feed_q)*feed_rate;
end
% Calculating the product composition
if feed_q >= 1
    for i = 1:length(x_eq)
       error(i) = abs(x_eq(i)-feed_Z);
       if error(i) < 0.001
          x_L = feed_Z;
          y_V = y_eq(i);
          break
       end
    end
elseif feed_q <= 0
    for i = 1:length(x_eq)
        error(i) = abs(y_eq(i)-feed_Z);
        if error(i) < 0.001
           y_V = feed_Z;
           x_V = x_eq(i);
           break
        end
    end
else
    for i = 1:length(x_eq)
        error(i) = abs((feed_q/(feed_q-1))*x_eq(i) - feed_Z/(feed_q-1) - y_eq(i));
        if error(i) < 0.001
            x_L = x_eq(i);
            y_V = y_eq(i);
            break
        end
    end
end

T0 = 300;                           % Initial guess (K)
if command1 ==1
    Bubble_T = @(T) TotalPressure(Water_Coeff(1),Water_Coeff(2),Water_Coeff(3),Methanol_Coeff(1),Methanol_Coeff(2),Methanol_Coeff(3) ... 
            ,T,feed_p_bar,x_L);     % Define Bubble Temp. function
    Dew_T = @(T) TotalPressureDew(Water_Coeff(1),Water_Coeff(2),Water_Coeff(3),Methanol_Coeff(1),Methanol_Coeff(2),Methanol_Coeff(3) ... 
            ,T,feed_p_bar,y_V);     % Define Dew Temp. function
    T_out_L = fsolve(Bubble_T,T0);  % Temperature of liquid product
    T_out_V = fsolve(Dew_T,T0);     % Temperature of vapor product
elseif command1 == 2
    Bubble_T = @(T) TotalPressure(Benzene_Coeff(1),Benzene_Coeff(2),Benzene_Coeff(3),Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3) ... 
            ,T,feed_p_bar,x_L);     % Define Bubble Temp. function
    Dew_T = @(T) TotalPressureDew(Benzene_Coeff(1),Benzene_Coeff(2),Benzene_Coeff(3),Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3) ... 
            ,T,feed_p_bar,y_V);     % Define Dew Temp. function
    T_out_L = fsolve(Bubble_T,T0);  % Temperature of liquid product (K)
    T_out_V = fsolve(Dew_T,T0);     % Temperature of vapor product (K)
elseif command1 == 3
    Bubble_T = @(T) TotalPressure(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),Octane_Coeff(1),Octane_Coeff(2),Octane_Coeff(3) ... 
            ,T,feed_p_bar,x_L);     % Define Bubble Temp. function
    Dew_T = @(T) TotalPressureDew(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),Octane_Coeff(1),Octane_Coeff(2),Octane_Coeff(3) ... 
            ,T,feed_p_bar,y_V);     % Define Dew Temp. function
    T_out_L = fsolve(Bubble_T,T0);  % Temperature of liquid product (K)
    T_out_V = fsolve(Dew_T,T0);     % Temperature of vapor product (K)
elseif command1 == 4
    Bubble_T = @(T) TotalPressure(Heptane_Coeff(1),Heptane_Coeff(2),Heptane_Coeff(3),CycloHexane_Coeff(1),CycloHexane_Coeff(2),CycloHexane_Coeff(3) ... 
            ,T,feed_p_bar,x_L);     % Define Bubble Temp. function
    Dew_T = @(T) TotalPressureDew(Heptane_Coeff(1),Heptane_Coeff(2),Heptane_Coeff(3),CycloHexane_Coeff(1),CycloHexane_Coeff(2),CycloHexane_Coeff(3) ... 
            ,T,feed_p_bar,y_V);     % Define Dew Temp. function
    T_out_L = fsolve(Bubble_T,T0);  % Temperature of liquid product (K)
    T_out_V = fsolve(Dew_T,T0);     % Temperature of vapor product (K)
elseif command1 == 5
    Bubble_T = @(T) TotalPressure(Propane_Coeff(1),Propane_Coeff(2),Propane_Coeff(3),Butane_Coeff(1),Butane_Coeff(2),Butane_Coeff(3) ... 
            ,T,feed_p_bar,x_L);     % Define Bubble Temp. function
    Dew_T = @(T) TotalPressureDew(Propane_Coeff(1),Propane_Coeff(2),Propane_Coeff(3),Butane_Coeff(1),Butane_Coeff(2),Butane_Coeff(3) ... 
            ,T,feed_p_bar,y_V);     % Define Dew Temp. function
    T_out_L = fsolve(Bubble_T,T0);  % Temperature of liquid product (K)
    T_out_V = fsolve(Dew_T,T0);     % Temperature of vapor product (K)
elseif command1 == 6
    Bubble_T = @(T) TotalPressure(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),CycloHexane_Coeff(1),CycloHexane_Coeff(2),CycloHexane_Coeff(3) ... 
            ,T,feed_p_bar,x_L);     % Define Bubble Temp. function
    Dew_T = @(T) TotalPressureDew(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),CycloHexane_Coeff(1),CycloHexane_Coeff(2),CycloHexane_Coeff(3) ... 
            ,T,feed_p_bar,y_V);     % Define Dew Temp. function
    T_out_L = fsolve(Bubble_T,T0);  % Temperature of liquid product (K)
    T_out_V = fsolve(Dew_T,T0);     % Temperature of vapor product (K)
end
clc;
%% Part6 - Calculating the preheater duty and steam mass
T_ref = 273.15;                 % Reference Temp. [=] K
T0 = 300;                       % Initial guess [K] 
if command1 ==1
        Bubble_T = @(T) TotalPressure(Water_Coeff(1),Water_Coeff(2),Water_Coeff(3),Methanol_Coeff(1),Methanol_Coeff(2),Methanol_Coeff(3) ... 
            ,T,feed_p_bar,feed_Z);              % Define Bubble Temp. function
        Dew_T = @(T) TotalPressureDew(Water_Coeff(1),Water_Coeff(2),Water_Coeff(3),Methanol_Coeff(1),Methanol_Coeff(2),Methanol_Coeff(3) ... 
            ,T,feed_p_bar,feed_Z);              % Define Dew Temp. function
        Temp_Bubble = fsolve(Bubble_T,T0);      % Determine Bubble Temperature (K)
        Temp_Dew = fsolve(Dew_T,T0);            % Determine Dew Temperature (K)
    if command3 == 1
      Cp_mix = feed_Z *Cp_Methanol_mol + (1-feed_Z)*Cp_Water_mol;       % Calculating Cp using mixture rule
      H1 = feed_rate * Cp_mix * (feed_T_K - T_ref);
      HL = feed_rate * Cp_mix * (Temp_Bubble - T_ref);
      HV = feed_rate * ((Temp_Dew - T_ref)*(Cp_mix) + (feed_Z*(hfg_Methanol*MW_Methanol)/1000 + (1-feed_Z)*(hfg_Water*MW_Water)/1000));
      if feed_q <= 0 
         H2 = HV; 
      elseif feed_q >= 1
          H2 = HL;
      else
          H2 = HV - feed_q*(HV - HL);
      end
    else
        Cp_mix = feed_Z *Cp_Methanol + (1-feed_Z)*Cp_Water;             % Calculating Cp using mixture rule
        H1 = feed_rate * Cp_mix * (feed_T_K - T_ref);
        HL = feed_rate * Cp_mix * (Temp_Bubble - T_ref);
        HV = feed_rate * ((Temp_Dew - T_ref)*(Cp_mix) + (feed_Z*hfg_Methanol + (1-feed_Z)*hfg_Water));
      if feed_q <= 0 
         H2 = HV; 
      elseif feed_q >= 1
          H2 = HL;
      else
          H2 = HV - feed_q*(HV - HL);
      end        
    end
    
elseif command1 == 2
        Bubble_T = @(T) TotalPressure(Benzene_Coeff(1),Benzene_Coeff(2),Benzene_Coeff(3),Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3) ... 
            ,T,feed_p_bar,feed_Z);              % Define Bubble Temp. function
        Dew_T = @(T) TotalPressureDew(Benzene_Coeff(1),Benzene_Coeff(2),Benzene_Coeff(3),Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3) ... 
            ,T,feed_p_bar,y_V);                 % Define Dew Temp. function
        Temp_Bubble = fsolve(Bubble_T,T0);      % Determine Bubble Temperature (K)
        Temp_Dew = fsolve(Dew_T,T0);            % Determine Dew Temperature (K)
    if command3 == 1
        Cp_mix = feed_Z * Cp_Toluene_mol + (1-feed_Z)*Cp_Benzene_mol;           % Calculating Cp using mixture rule
        H1 = feed_rate * Cp_mix * (feed_T_K - T_ref);
        HL = feed_rate * Cp_mix * (Temp_Bubble - T_ref);
        HV = feed_rate * ((Temp_Dew - T_ref)*(Cp_mix) + (feed_Z*(hfg_Toluene*MW_Toluene)/1000 + (1-feed_Z)*(hfg_Benzene*MW_Benzene)/1000));
      if feed_q <= 0 
         H2 = HV; 
      elseif feed_q >= 1
          H2 = HL;
      else
          H2 = HV - feed_q*(HV - HL);
      end        
    else
        Cp_mix = feed_Z * Cp_Toluene + (1-feed_Z)*Cp_Benzene;                   % Calculating Cp using mixture rule
        H1 = feed_rate * Cp_mix * (feed_T_K - T_ref);
        HL = feed_rate * Cp_mix * (Temp_Bubble - T_ref);
        HV = feed_rate * ((Temp_Dew - T_ref)*(Cp_mix) + (feed_Z*hfg_Toluene + (1-feed_Z)*hfg_Benzene));
      if feed_q <= 0 
         H2 = HV; 
      elseif feed_q >= 1
          H2 = HL;
      else
          H2 = HV - feed_q*(HV - HL);
      end        
    end
    
elseif command1 == 3
        Bubble_T = @(T) TotalPressure(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),Octane_Coeff(1),Octane_Coeff(2),Octane_Coeff(3) ... 
            ,T,feed_p_bar,feed_Z);              % Define Bubble Temp. function
        Dew_T = @(T) TotalPressureDew(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),Octane_Coeff(1),Octane_Coeff(2),Octane_Coeff(3) ... 
            ,T,feed_p_bar,y_V);                 % Define Dew Temp. function
        Temp_Bubble = fsolve(Bubble_T,T0);      % Determine Bubble Temperature (K)
        Temp_Dew = fsolve(Dew_T,T0);            % Determine Dew Temperature (K)
    if command3 == 1
        Cp_mix = feed_Z * Cp_Toluene_mol + (1-feed_Z)*Cp_Octane_mol;                    % Calculating Cp using mixture rule
        H1 = feed_rate * Cp_mix * (feed_T_K - T_ref);
        HL = feed_rate * Cp_mix * (Temp_Bubble - T_ref);
        HV = feed_rate * ((Temp_Dew - T_ref)*(Cp_mix) + (feed_Z*(hfg_Toluene*MW_Toluene)/1000 + (1-feed_Z)*(hfg_Octane*MW_Octane)/1000));
      if feed_q <= 0 
         H2 = HV; 
      elseif feed_q >= 1
          H2 = HL;
      else
          H2 = HV - feed_q*(HV - HL);
      end        
    else
        Cp_mix = feed_Z * Cp_Toluene + (1-feed_Z)*Cp_Octane;                            % Calculating Cp using mixture rule
        H1 = feed_rate * Cp_mix * (feed_T_K - T_ref);
        HL = feed_rate * Cp_mix * (Temp_Bubble - T_ref);
        HV = feed_rate * ((Temp_Dew - T_ref)*(Cp_mix) + (feed_Z*hfg_Toluene + (1-feed_Z)*hfg_Octane));
      if feed_q <= 0 
         H2 = HV; 
      elseif feed_q >= 1
          H2 = HL;
      else
          H2 = HV - feed_q*(HV - HL);
      end        
    end
    
elseif command1 == 4
        Bubble_T = @(T) TotalPressure(Heptane_Coeff(1),Heptane_Coeff(2),Heptane_Coeff(3),CycloHexane_Coeff(1),CycloHexane_Coeff(2),CycloHexane_Coeff(3) ... 
            ,T,feed_p_bar,feed_Z);                  % Define Bubble Temp. function
        Dew_T = @(T) TotalPressureDew(Heptane_Coeff(1),Heptane_Coeff(2),Heptane_Coeff(3),CycloHexane_Coeff(1),CycloHexane_Coeff(2),CycloHexane_Coeff(3) ... 
            ,T,feed_p_bar,y_V);                     % Define Dew Temp. function
        Temp_Bubble = fsolve(Bubble_T,T0);          % Determine Bubble Temperature (K)
        Temp_Dew = fsolve(Dew_T,T0);                % Determine Dew Temperature (K)
    if command3 == 1
        Cp_mix = feed_Z * Cp_CycloHexane_mol + (1-feed_Z)*Cp_Heptane_mol;                   % Calculating Cp using mixture rule
        H1 = feed_rate * Cp_mix * (feed_T_K - T_ref);
        HL = feed_rate * Cp_mix * (Temp_Bubble - T_ref);
        HV = feed_rate * ((Temp_Dew - T_ref)*(Cp_mix) + (feed_Z*(hfg_CycloHexane*MW_CycloHexane)/1000 + (1-feed_Z)*(hfg_Heptane*MW_Heptane)/1000));
      if feed_q <= 0 
         H2 = HV; 
      elseif feed_q >= 1
          H2 = HL;
      else
          H2 = HV - feed_q*(HV - HL);
      end        
    else
        Cp_mix = feed_Z * Cp_CycloHexane + (1-feed_Z)*Cp_Heptane;                           % Calculating Cp using mixture rule
        H1 = feed_rate * Cp_mix * (feed_T_K - T_ref);
        HL = feed_rate * Cp_mix * (Temp_Bubble - T_ref);
        HV = feed_rate * ((Temp_Dew - T_ref)*(Cp_mix) + (feed_Z*hfg_CycloHexane + (1-feed_Z)*hfg_Heptane));
      if feed_q <= 0 
         H2 = HV; 
      elseif feed_q >= 1
          H2 = HL;
      else
          H2 = HV - feed_q*(HV - HL);
      end        
    end
    
elseif command1 == 5
        Bubble_T = @(T) TotalPressure(Propane_Coeff(1),Propane_Coeff(2),Propane_Coeff(3),Butane_Coeff(1),Butane_Coeff(2),Butane_Coeff(3) ... 
            ,T,feed_p_bar,feed_Z);                  % Define Bubble Temp. function
        Dew_T = @(T) TotalPressureDew(Propane_Coeff(1),Propane_Coeff(2),Propane_Coeff(3),Butane_Coeff(1),Butane_Coeff(2),Butane_Coeff(3) ... 
            ,T,feed_p_bar,y_V);                     % Define Dew Temp. function
        Temp_Bubble = fsolve(Bubble_T,T0);          % Determine Bubble Temperature (K)
        Temp_Dew = fsolve(Dew_T,T0);                % Determine Dew Temperature (K)
    if command3 == 1
        Cp_mix = feed_Z * Cp_Propane_mol + (1-feed_Z)*Cp_Butane_mol;                        % Calculating Cp using mixture rule
        H1 = feed_rate * Cp_mix * (feed_T_K - T_ref);
        HL = feed_rate * Cp_mix * (Temp_Bubble - T_ref);
        HV = feed_rate * ((Temp_Dew - T_ref)*(Cp_mix) + (feed_Z*(hfg_Propane*MW_Propane)/1000 + (1-feed_Z)*(hfg_Butane*MW_Butane)/1000));
      if feed_q <= 0 
         H2 = HV; 
      elseif feed_q >= 1
          H2 = HL;
      else
          H2 = HV - feed_q*(HV - HL);
      end        
    else
        Cp_mix = feed_Z * Cp_Propane + (1-feed_Z)*Cp_Butane;                                % Calculating Cp using mixture rule
        H1 = feed_rate * Cp_mix * (feed_T_K - T_ref);
        HL = feed_rate * Cp_mix * (Temp_Bubble - T_ref);
        HV = feed_rate * ((Temp_Dew - T_ref)*(Cp_mix) + (feed_Z*hfg_Propane + (1-feed_Z)*hfg_Butane));
      if feed_q <= 0 
         H2 = HV; 
      elseif feed_q >= 1
          H2 = HL;
      else
          H2 = HV - feed_q*(HV - HL);
      end        
    end
    
elseif command1 == 6
        Bubble_T = @(T) TotalPressure(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),CycloHexane_Coeff(1),CycloHexane_Coeff(2),CycloHexane_Coeff(3) ... 
            ,T,feed_p_bar,feed_Z);                  % Define Bubble Temp. function
        Dew_T = @(T) TotalPressureDew(Toluene_Coeff(1),Toluene_Coeff(2),Toluene_Coeff(3),CycloHexane_Coeff(1),CycloHexane_Coeff(2),CycloHexane_Coeff(3) ... 
            ,T,feed_p_bar,y_V);                     % Define Dew Temp. function
        Temp_Bubble = fsolve(Bubble_T,T0);          % Determine Bubble Temperature (K)
        Temp_Dew = fsolve(Dew_T,T0);                % Determine Dew Temperature (K)
    if command3 == 1
        Cp_mix = feed_Z * Cp_Toluene_mol + (1-feed_Z)*Cp_CycloHexane_mol;                       % Calculating Cp using mixture rule
        H1 = feed_rate * Cp_mix * (feed_T_K - T_ref);
        HL = feed_rate * Cp_mix * (Temp_Bubble - T_ref);
        HV = feed_rate * ((Temp_Dew - T_ref)*(Cp_mix) + (feed_Z*(hfg_Toluene*MW_Toluene)/1000 + (1-feed_Z)*(hfg_CycloHexane)/1000));
      if feed_q <= 0 
         H2 = HV; 
      elseif feed_q >= 1
          H2 = HL;
      else
          H2 = HV - feed_q*(HV - HL);
      end  
      
    else
       Cp_mix = feed_Z * Cp_Toluene + (1-feed_Z)*Cp_CycloHexane;                                % Calculating Cp using mixture rule
       H1 = feed_rate * Cp_mix * (feed_T_K - T_ref);
       HL = feed_rate * Cp_mix * (Temp_Bubble - T_ref);
       HV = feed_rate * ((Temp_Dew - T_ref)*(Cp_mix) + (feed_Z*hfg_Toluene + (1-feed_Z)*hfg_CycloHexane));
      if feed_q <= 0 
         H2 = HV; 
      elseif feed_q >= 1
          H2 = HL;
      else
          H2 = HV - feed_q*(HV - HL);
      end
    end
end
clc;

%% Part7 - Showing results 
% Result of part 5
Liquid_Product_abs_Temperature = T_out_L';
Vapor_Product_abs_Temperature = T_out_V';
Liquid_Amount = L_out';
Vapor_Amount = V_out';
Liquid_Composition = x_L';
Vapor_Composition = y_V';
disp('The output properties of column is:(Temperatures are in Kelvin) ')
table1 = table(Liquid_Amount,Vapor_Amount,Liquid_Composition,Vapor_Composition,Liquid_Product_abs_Temperature,Vapor_Product_abs_Temperature);
disp(table1)
% Result of part6
delta_H = H2 - H1;
m_vapor = delta_H / hfg;
disp('The preheater duty is (J): '); disp(delta_H)
disp('The rate of steam is (kg/hr): '); disp(m_vapor)
% Result of part3
figure
subplot(2,1,1)
plot(x_eq,y_eq);hold on
plot([0,1],[0,1])
plot(x_L,y_V,'B.','MarkerSize',20)
plot([feed_Z,x_L],[feed_Z,y_V])
hold off; grid on; grid minor
xlabel('x_i'); ylabel('y_i')
subplot(2,1,2)
plot(x_eq,T_Bub); hold on 
plot(y_eq,T_Bub)
hold off; grid on; grid minor
xlabel('x_i , y_i'); ylabel('T (Celcius)')
% Bonus
x_L_vector = linspace(0,1,11);
y_V_vector = linspace(0,1,11);
q_Bonus = linspace(0,1,11);
for i = 1:length(q_Bonus)
    for j = 1:length(x_eq)
        error(i) = abs((q_Bonus(i)/(q_Bonus(i)-1))*x_eq(j) - feed_Z/(q_Bonus(i)-1) - y_eq(j));
        if error(i) < 0.001
            x_L_vector(i) = x_eq(j);
            y_V_vector(i) = y_eq(j);
            break
        end
    end
end
figure
plot(x_eq,y_eq); hold on 
plot([0,1],[0,1])
plot(x_L_vector(1),y_V_vector(1),'B.','MarkerSize',20);   plot([feed_Z,x_L_vector(1)],[feed_Z,y_V_vector(1)])
plot(x_L_vector(2),y_V_vector(2),'B.','MarkerSize',20);   plot([feed_Z,x_L_vector(2)],[feed_Z,y_V_vector(2)])
plot(x_L_vector(3),y_V_vector(3),'B.','MarkerSize',20);   plot([feed_Z,x_L_vector(3)],[feed_Z,y_V_vector(3)])
plot(x_L_vector(4),y_V_vector(4),'B.','MarkerSize',20);   plot([feed_Z,x_L_vector(4)],[feed_Z,y_V_vector(4)])
plot(x_L_vector(5),y_V_vector(5),'B.','MarkerSize',20);   plot([feed_Z,x_L_vector(5)],[feed_Z,y_V_vector(5)])
plot(x_L_vector(6),y_V_vector(6),'B.','MarkerSize',20);   plot([feed_Z,x_L_vector(6)],[feed_Z,y_V_vector(6)])
plot(x_L_vector(7),y_V_vector(7),'B.','MarkerSize',20);   plot([feed_Z,x_L_vector(7)],[feed_Z,y_V_vector(7)])
plot(x_L_vector(8),y_V_vector(8),'B.','MarkerSize',20);   plot([feed_Z,x_L_vector(8)],[feed_Z,y_V_vector(8)])
plot(x_L_vector(9),y_V_vector(9),'B.','MarkerSize',20);   plot([feed_Z,x_L_vector(9)],[feed_Z,y_V_vector(9)])
plot(x_L_vector(10),y_V_vector(10),'B.','MarkerSize',20); plot([feed_Z,x_L_vector(10)],[feed_Z,y_V_vector(10)])
plot(x_L_vector(11),y_V_vector(11),'B.','MarkerSize',20); plot([feed_Z,x_L_vector(11)],[feed_Z,y_V_vector(11)])
hold off; grid on; grid minor
xlabel('x_i'); ylabel('y_i'); title('x-y diagram and differents q from 0 to 1')