%% ----------------------- Import Flight data  ------------------------

FD = open('FTISxprt-20190308_125059.mat') ;     
Time = FD.flightdata.time.data ;    %Time in (mili)seconds

%% --------------------------- Constants ------------------------------

g = 9.81 ;          % [m/s^2] gravitational constant
gamma = 1.4 ;       % spec. heat ratio
R = 287.05 ;        % [J/kg*K] - spec. gas constant
S = 30.00 ;	        % [m^2] wing area
c = 2.0569;         % mean aerodynamic cord [m]
E_d = 0.686 ;       % [m] characteristic diameter JT15D-4 engine
bpr = 2.6 ;         % by-pass ratio
mf = 35.245 ;       % [kg/s] - total mass flow
mff = 0.177;        % [kg/s] - engine fuel flow
lambda = -0.0065;   % [K/m] - temperature gradient in ISA
Ws = 60500 ;        % [N] - standard aircraft weight
Cmde = -1.1642 ;    % Elevator effectiveness
CmTc = -0.0064 ;    % dimensionless thrust moment arm

%% -------------------------- datasheet ------------------------------

hp0 = [11900, 12050, 12090, 12270, 11320, 10700, 10250]' ;  % [ft]
Vc = [160, 149, 138, 130, 170, 181, 190]' ;                 % [kts]
TAT = [-11.8, -12.8, -13.8, -14.2, -10, -8.1, -6.3]' ;      % [C]
mfl = [391, 388, 388, 385, 402, 411, 418]' ;                % [lbs/hr]
mfr = [422, 418, 419, 415, 434, 445, 451]' ;                % [lbs/hr]
mffs = [0.048, 0.048, 0.048, 0.048, 0.048, 0.048, 0.048]' ; % [kg/s]
Wfu = [764, 779, 804, 819, 840, 856, 871]' ;                % [lbs]
de = [0.2, -0.2, -0.7, -1.1, 0.5, 0.85, 1.1]' ;             % [deg]
Fe = [1, -16, -30, -36, 26, 55, 80]' ;                      % [N]
alpha = [5.2, 6.2, 7.2, 8.2, 4.4, 3.8, 3.4]' ;              % [deg]

% sea-level ISA parameters
[T_isa, a_isa, P_isa, rho_isa] = atmoscoesa(0) ;

%% ---------------------- Conversion of units ------------------------

% measured total temp.
TAT = TAT + 273.15 ;    % [K]

% pressure altitude
hp0 = hp0./3.2808399 ;  % [m]

% fuel flow left engine
mfl = mfl.*(0.45359237/3600) ;   % [kg/s]

% fuel flow right engine
mfr = mfr.*(0.45359237/3600) ;   % [kg/s]

% calibrated airspeed 
Vc = Vc.*0.51444 ; %[m/s]

%% --------------- Reduction of the measured airspeed -----------------

[Mach, Ve, rho] = speedred(Vc, TAT, hp0, T_isa, P_isa, rho_isa, g, R, gamma, lambda) ;

%% --------------------- Thrust calculations --------------------------

[T, Ts] = thrustcalc(TAT, hp0, lambda, Mach, mfl, mfr, mffs) ;

%% ---------- Reduction of the non-standard aircraft mass -------------

[Mass, W, Ve_r, Mass_ramp] = massred(Wfu, Ws, Ve, g) ;

%% ---------- Reduction of the non-standard engine thrust -------------

de_r = elevred(T, Ts, de, CmTc, Cmde, rho, rho_isa, Vc, E_d) ;

%% ------------ Drawing the reduced elevator trim curve ---------------

trimcurve(Ve_r, de_r)

%% -------------- Reduced elevator control force curve ----------------

Fe_r = controlcurve(Fe, Ws, W, Ve_r) ;


%% --------------------- Elevator effectiveness -----------------------

Cmde_m = eleveff(Mass_ramp, rho, S, c, g) ;


%% --------------------- Longitudinal stability -----------------------

Cma_m = longstab(alpha, de, Cmde_m) ;

%% ---------------------- Conversion functions ------------------------

function pounds = KG2P(kg)

    pounds = kg / 0.453592;
    
end

function kg = P2KG(pounds)

    kg = pounds * 0.453592;
    
end

function [Mach, Ve, rho] = speedred(Vc, TAT, hp0, T_isa, P_isa, rho_isa, g, R, gamma, lambda)

    for i=1:7

        V_c = Vc(i) ;           % [m/s] - calibrated airspeed
        T_m = TAT(i) ;          % [K] - measured air temp.
        h_p = hp0(i);           % [m] - ISA pressure altitude

        % static pressure
        P = P_isa*(1 + ((lambda*h_p)/T_isa))^(-g/(lambda*R)) ;

        % Mach number
        Mch = sqrt((2/(gamma - 1))*((1 + (P_isa/P)*((1 + ((gamma - 1)*rho_isa*(V_c)^2)/(2*gamma*P_isa))^(gamma/(gamma - 1)) - 1))^((gamma - 1)/gamma) - 1)) ;

        M(i) = Mch ;

        % correct for ram rise
        T = T_m/(1+((gamma - 1)/2)*Mch^2) ;

        % density based on ideal gas law
        rho = P/(R*T) ;

        % local speed of sound
        a = sqrt(gamma*R*T) ;

        % true airspeed
        V_t = Mch*a ;

        % equivalent airspeed
        V_e = V_t*sqrt(rho/rho_isa) ;

        Ve(i) = V_e ;

    end

    Mach = M' ;
    Ve = Ve' ;
    
end

function [T, Ts] = thrustcalc(TAT, hp0, lambda, Mach, mfl, mfr, mffs)

    % ISA Temperature
    T_ISA = 288.15 + (lambda*hp0) ; % [K]

    % Temp. difference
    dT = TAT - T_ISA ;              % [K]

    % Merging data set - correct format
    data = [hp0 Mach dT mfl mfr] ;
    dlmwrite('matlab.dat',data,'delimiter',' ') 

    % Execute thrust calculations
    !Thrust.exe & ;

    pause(1) ;

    T1 = load('thrust.dat') ;  % load thrust force per engine

    T = sum(T1,2) ;

    % STANDARDIZED VERSION

    datas = [hp0 Mach dT mffs mffs] ;
    dlmwrite('matlab.dat',datas,'delimiter',' ') 

    % Execute thrust calculations
    !Thrust.exe & ;

    pause(1) ;

    T2 = load('thrust.dat') ;   % load thrust force per engine

    Ts = sum(T2,2) ;
    
end

function [Mass, W, Ve_r, Mass_ramp] = massred(Wfu, Ws, Ve, g)

    %START MASS IN POUNDS
    BEM  =  9165 ;    %Provided Basic Empty Mass in POUNDS
    FUEL0=  3850 ;    %Total Fuel on T=0

    Mass_pax  = KG2P([82,92,56,63,75,75,76,77,75]) ; 
    Mass_bag  = KG2P([0,0,0]) ;                       
    Mass_PAY  = sum(Mass_pax) + sum(Mass_bag) ;   %Total payload in POUNDS

    Mass_ramp =  BEM + FUEL0 + Mass_PAY ;

    for j=1:7

        % aircraft mass
        Mass_t = Mass_ramp - Wfu(j) ; % [lbs]

        Mass_t = P2KG(Mass_t) ; % [kg]

        Mass(j) = Mass_t ;

        % aircraft weight
        Wght = Mass_t*g ;

        W(j) = Wght;    % [N]

        % reduced equivalent airspeed
        Ver = Ve(j)*sqrt(Ws/Wght) ; %[m/s]
        Ve_r(j) = Ver ;

    end

    Mass = Mass' ;
    W = W' ;
    Ve_r = Ve_r' ;
    
end

function de_r = elevred(T, Ts, de, CmTc, Cmde, rho, rho_isa, Vc, E_d)

    A = ((E_d/2)^(2))*pi ; % area of the engine inlet

    Tc = T./(0.5*rho*(Vc.^(2))*A) ;
    Tcs = Ts./(0.5*rho_isa*(Vc.^(2))*A) ;

    for b=1:7

        % reduced equivalent airspeed
        der = de(b) - ((CmTc/Cmde)*(Tcs(b) - Tc(b))); 

        de_r(b) = der ;

    end

    de_r = de_r' ;

end

function trimcurve(Ve_r, de_r)

    D = sort(Ve_r) ;
    E = sort(de_r) ;

    figure(1)
    plot(D,E,'.b')
    p = polyfit(D,E,2) ;
    f = polyval(p,D);
    hold on
    plot(D,f,'-r')
    ax = gca ;
    ax.YDir = 'reverse' ;
    plot(65:95, zeros(31, 1),'k')

    legend('data points','Elevator trim curve')
    text(70, -1.25, 'Configuration: (Gear up, flaps up)');
    text(75, -0.5, 'Cm_{\alpha} < 0');

    title('Elevator trim curve of the Cessna Citation - II PH-LAB')
    xlabel('Reduced equivalent airspeed V_{e}')
    ylabel('Reduced elevator deflection \delta_{e}')

    grid on

end

function Fe_r = controlcurve(Fe, Ws, W, Ve_r)

    Fe_aer = Fe ;

    for s=1:7

        Fe_aer_temp = Fe_aer(s)*sqrt(Ws/W(s)) ;

        Fe_aer_r(s) = Fe_aer_temp ;

    end

    Fe_r = Fe_aer_r' ;

    C = sort(Ve_r) ;
    B = sort(Fe_r) ;

    figure(2)
    plot(C,B,'.b')
    p2 = polyfit(C,B,2) ;
    f2 = polyval(p2,C);
    hold on
    plot(C,f2,'-r')
    ax = gca ;
    ax.YDir = 'reverse' ;
    plot(65:95, zeros(31, 1),'k')

    legend('data points','Elevator control force curve')
    text(70, 40, 'Configuration: (Gear up, flaps up)');

    title('Elevator control force curve of the Cessna Citation - II PH-LAB')
    xlabel('Reduced equivalent airspeed V_{e}')
    ylabel('Reduced elevator control force F_{e}')

    grid on

end

function Cmde_m = eleveff(Mass_ramp, rho, S, c, g)

    Mxcg = Mass_ramp - 905 ;    % [lbs]
    Mxcg = P2KG(Mxcg) ;         % [kg]
    Wxcg = Mxcg*g ;             % [N]
    VIAS = 160*0.51444 ;        % [m/s]

    dxcg = (282.01 - 280.04)*0.0254 ;
    de_dxcg = (0.2 - (-0.3))*pi/180 ;
    C_N = Wxcg/(0.5*rho*(VIAS^(2))*S) ;

    Cmde_m = -(1/de_dxcg)*C_N*(dxcg/c) ;
    
end

function Cma_m = longstab(alpha, de, Cmde_m)

    A = alpha ;
    F = de ;

    figure(3)
    plot(A,F,'.b')
    p3 = polyfit(A,F,1) ;
    f3 = polyval(p3,A);
    hold on
    plot(A,f3,'-r')
    ax = gca ;
    ax.YDir = 'reverse' ;
    plot(3:8, zeros(6, 1),'k')

    legend('data points', 'Elevator trim curve', 'Location','northwest')
    text(70, 40, 'Configuration: (Gear up, flaps up)');

    title('Elevator trim curve w.r.t. AoA of the Cessna Citation - II PH-LAB')
    xlabel('Angle of Attack')
    ylabel('Elevator deflection \delta_{e}')

    grid on
    
    grad = (f3(7)-f3(1))/(A(7)-A(1)) ;
    
    Cma_m = -grad*Cmde_m ;

end
