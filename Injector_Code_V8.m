clc; clear;

% OBJECTIVES:
%   Spray Angle = 45 deg
%   Minimize:
%       Sauter Mean Diameter (SMD)
%       Vaporization Distance (VD)
%       Transit Time (t)
%       Combustion Chamber Length
%       Combustion Chamber Diameter
%
% UNITS: lbf, lbm, in, °F, °R, Btu, psi, s
%_____________________________________________________________________________

% VARIABLE INPUTS
    %Assumed pressure drop of 20%
    %Assumed values
    p_drop_f = 60;      % Pressure drop of fuel [psi]
    p_drop_ox = 60;     % Pressure drop of oxidizer [psi]
    Cd_f = .70;         % Discharge coefficient of fuel orifice (assumed from chart)
    Cd_ox = .65;         % Discharge coefficient of oxidizer orifice (assumed from chart)
    conc_f = 1;         % Ethanol-water Concentration

% CONSTANTS
    rho_f = .0285*conc_f + 0.036*(1-conc_f);          % Density of fuel (% eth + % water) [lb/in^3]    
    rho_ox = .0412;                                   % Density of oxidizer [lb/in^3] (=1141 %kg/m^3)   
    mu_f = 22.43E10-6;         % fuel dynamic viscosity @ STD [lbf*s/ft^2]
    R = 2400;                  % Gas constant of combustion chamber gas [in-lb/lb°R](Assumed)
    sigma_f = 0.0003237242;    % surface tension of fuel [psi]
    T_inj_ox = 97.848;         % injection temperature of LoX [°R]
    k_ox = 4.95E-5;            % thermal conductivity of LoX [Btu-in/(s-in^2-°R)]
    k_f = 3.1845E-5;           % thermal conductivity of methane @ 100K [Btu-in/(s-ft^2-°R)]
    k_c = 4.95E-5;             %Thermal Conductivity of chamber gas @ X temp [Btu-in/(s-ft^2-°R)]
    Cp_ox = 2.64;              % gas specific heat [Btu/lb°R]
    Cp_c = 2.64;               % combustion chamber gas specific heat
    gam = 1.2;                 % specific heat ratio (assumed)
    Pr_c = (4*gam)/((9*gam)-5);% combustion chamber gas prandtle number(approximated)

% PINTLE DIMENSIONS
    %A_ox = .012;       % Exit (orifice) area of oxidizer [in^2]
    %A_f = 0.5;         % Exit (orifice) area of fuel [in^2]
    Sl_t = .125;        % Pintle sleeve thickness [in]
    L_open = .125;        % pintle opening distance [in]
    theta_pt = 0;       % pintle tip angle [deg]
    ann_t = .07;        % annular gap thickness [in]
    %theta_UTL = unitless angle(in radians???)
    D_pt = 1;           % Pintle opening distance diameter [in]
    r_post = D_pt;      % length pintle protudes out [in]
%     L_min = lower limit of pintle opening distance [in]
%     A_min = minimum sleeve area between the sleeve and pintle tip [in]

% INJECTOR DIMENSIONS
    A_pt = L_open*pi*(D_pt - (2*Sl_t))
%     L_open = (L_min/cos(theta_pt));
%     A_min = pi*(D_pt - (2*Sl_t) - (L_min*sin(theta_pt)));
    A_ann = pi*((ann_t + (D_pt/2))^2 - (D_pt/2)^2)

% CHAMBER DIMENSIONS(ASSUMED)
     A_ch = 5.75;        % chamber cross sectional area [in^2]
     A_th =  1;          % throat area [in^2]
%_____________________________________________________________________________

% FLUID VELOCITY
    % dependent on pressure, TBD
    U_f = Cd_f * sqrt((2*p_drop_f)/rho_f);                  % Velocity of fuel [in/s]
    U_ox = Cd_ox * sqrt((2*p_drop_ox)/rho_ox);   % Velocity of oxidizer [in/s]

% MASS FLOW RATE
    m_dot_ox = (rho_ox) *(A_pt) * U_ox ;                 % oxidizer mass flow rate [lb/s]
    m_dot_f = (rho_f) *(A_ann) * U_f;                      % fuel mass flow rate [lb/s]
    OF = m_dot_ox/m_dot_f;                               % Oxidizer-Fuel Ratio
    m_dot_ch = 20;%(rho_f) *(A_ch) * U_f ;                   % chamber mass flow rate [lb/s]

% TOTAL MOMENTUM RATIO (TMR)
    TMR = ((m_dot_ox*U_ox)*cos(theta_pt))/((m_dot_f*U_f)+(m_dot_ox*U_ox)*sin(theta_pt));

%Spray angle 
    Spray_Angle = acos(1/(1 + TMR)) * (180/pi);          % [rad to degrees]
%_____________________________________________________________________________

% SAUTER MEAN DIAMETER (SMD)

    We = (rho_ox * L_open*(U_ox - U_f)^2)/sigma_f;           % Weber number
    theta_UTL = (90 - theta_pt)/90;                          % Non-Dimensional pintle tip angle
    q = 3.455 - (0.225*theta_UTL);
    SMD = L_open * (theta_UTL^-1) * exp(4.0 - q*(We^.1));    % SMD [in]
%_____________________________________________________________________________

% VAPORIZATION CALCULATIONS
        T_0 = T_inj_ox + U_ox^2/(2*Cp_ox);       % oxidizer stagnation temperature [°R]
        T_f = 567;                               % fuel stagnation temperature (assumed) [°R]
        h_fg = 364;                              % fuel heat of vaporization (assumed) [Btu/lb]
        Tc = 5400;                               % combustion chamber gas temperature
  % Drag (frictional drag between drop and gas flow)
        B = Cp_c*(T_f-T_0)/h_fg;
        Drag = (9/2)*(Pr_c/B); %Drag Coefficient
    % Vaporization Distance (VD)(in)
        VDC =  ((gam+1)/2)^((gam+1)/2*(gam-1));  %variable for droplet vaportization
        VD = ((SMD/2)^2)*((((U_f+U_ox)/2)/sqrt(gam*R*Tc)) + ((3/VDC)*(A_th/A_ch)*(Drag/10))*((Cp_ox*rho_ox)/k_c)*((sqrt(gam*R*Tc)/log(1+B))*(1/(2*Drag))));
    % Required Chamber Length(in)
        L_ch = VD*cosd(Spray_Angle/2) + L_open + r_post;
    %Required Chamber Diameter(in)
        D_ch = 2*VD*sind(Spray_Angle/2) + D_pt;
    % Transit Time (the time that the drop resides in  combustion chamber)
    % (s)
        t = L_ch/sqrt(gam*R*T_0)*((gam+1)/2)^((gam+1)/(2*(gam-1)))*(A_ch/A_th); 
%_____________________________________________________________________________

%Thrust Calculation(Rough)
%T = m_dot_tot * Ve
%Sum fuel and oxidizer m dots to obtain (lb/s)
%Ve = exhaust velocity(ft/s)

Z = [U_f; U_ox; m_dot_f; m_dot_ox; TMR; Spray_Angle; SMD; VD; Drag; L_ch;D_ch; t; OF];

%Matrix Outputs
    Inj_Dimensions = [A_pt, A_ann,Sl_t,L_open, D_pt, r_post, theta_pt]
    Inj_Perf = [U_f,U_ox,m_dot_f, m_dot_ox Spray_Angle, SMD, VD, Drag, t]
    CC_Dimenstions = [L_ch, D_ch]

%Table of performance and combustion chamber outputs
T = table(U_f', U_ox', m_dot_f', m_dot_ox', TMR', Spray_Angle',SMD',VD',Drag',L_ch',D_ch',t',OF','VariableNames',{'Fuel Velocity(in/s)','Oxidizer Velocty(in/s)','Fuel mdot(lb/s)','Oxidizer mdot(lb/s)','TMR','Spray Angle(deg)','SMD(in)','Vaporization Dist (in)','Drag','Req. Chamber Length (in)','Chamber Diameter(in)','Transit Time (s)','OF Ratio'})

