g = 9.81 ; % gravitational constant
gamma = 1.4 ;   % spec. heat ratio
R = 287.05 ; % [J/kg*K] - spec. gas constant
bpr = 2.6 ;     % by-pass ratio
mf = 35.245 ;    % [kg/s] - total mass flow

% -------------------------------------------------------------------

% datasheet
hp0 = transpose([6990, 7000, 8010, 8000, 8010, 8980]) ;   % [ft]
Vc = transpose([241, 219, 192, 161, 137, 119]) ;         % [kts]
TAT = transpose([1.2, -1.2, -2.2, -6.5, -6.2, -7.0]) ;  % [C]
mfl = transpose([720, 603, 492, 464, 387, 373]) ;       % [lbs/hr]
mfr = transpose([756, 643, 533, 509, 426, 409]) ;       % [lbs/hr]

% ---------------------------------------------------------------------

% measured total temp.
TAT = TAT + 273.15 ;    % [K]

% pressure altitude
hp0 = hp0./3.2808399 ;      % [m]

% fuel flow left engine
mfl = mfl*(0.45359237/3600) ;   % [kg/s]

% fuel flow right engine
mfr = mfr*(0.45359237/3600) ;   % [kg/s]

% Mach number
for i=1:6
    
    m_f1 = mfl(i) ;         % [kg/s] - fuel flow left engine
    m_f2 = mfr(i) ;         % [kg/s] - fuel flow right engine
    V_c = Vc(i)*0.51444 ;    % [m/s] - calibrated airspeed
    T_m = TAT(i) ;          % [K] - measured air temp.

    h_p = hp0(i);             % [m] - ISA pressure altitude
    [T_isa, a_isa, P_isa, rho_isa] = atmoscoesa(h_p) ;  % ISA parameters

    % ---------- Reduction of the airspeed -------------

    % static pressure
    P = P_isa*(1 + ((bpr*h_p)/T_isa))^(g/(bpr*R)) ;

    % Mach number
    Mch = sqrt((2/(gamma - 1))*((1 + (P_isa/P)*((1 + ((gamma - 1)*rho_isa*(V_c)^2)/(2*gamma*P_isa))^(gamma/(gamma - 1)) - 1))^((gamma - 1)/gamma) - 1)) ;

    M(i) = Mch;

end

Mach = transpose(M) ;

% ISA Temperature
T_ISA = 288.15 + ((-0.0065)*hp0) ;      % [K]

% Temp. difference
dT = TAT - T_ISA ;                      % [K]

% ---------------------------------------------------------------------

% Merging data set - correct format
data = [hp0 Mach dT mfl mfr] ;
data = round(data*1000)/1000 ;
dlmwrite('matlab.dat',data,'delimiter',' ') 
load matlab.dat   % load the file
