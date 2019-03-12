g = 9.81 ;          % gravitational constant
gamma = 1.4 ;       % spec. heat ratio
R = 287.05 ;        % [J/kg*K] - spec. gas constant
bpr = 2.6 ;         % by-pass ratio
mf = 35.245 ;       % [kg/s] - total mass flow
lambda = -0.0065;   % temperature gradient in ISA [K/m]

% -------------------------------------------------------------------

% datasheet
hp0 = [6990, 7000, 8010, 8000, 8010, 8980]' ; % [ft]
Vc = [241, 219, 192, 161, 137, 119]' ;        % [kts]
TAT = [1.2, -1.2, -2.2, -6.5, -6.2, -7.0]' ;  % [C]
mfl = [720, 603, 492, 464, 387, 373]' ;       % [lbs/hr]
mfr = [756, 643, 533, 509, 426, 409]' ;       % [lbs/hr]

[T_isa, a_isa, P_isa, rho_isa] = atmoscoesa(0) ;      % ISA parameters

% --------------------------------------------------

% measured total temp.
TAT = TAT + 273.15 ;    % [K]

% pressure altitude
hp0 = hp0/3.2808399 ;  % [m]

% fuel flow left engine
mfl = mfl*(0.45359237/3600) ;   % [kg/s]

% fuel flow right engine
mfr = mfr*(0.45359237/3600) ;   % [kg/s]

% ---------- Reduction of the airspeed -------------

for i=1:6
    
    V_c = Vc(i)*0.51444 ;   % [m/s] - calibrated airspeed
    T_m = TAT(i) ;          % [K] - measured air temp.
    h_p = hp0(i);           % [m] - ISA pressure altitude

    % static pressure
    P = P_isa*(1 + ((lambda*h_p)/T_isa))^(g/(lambda*R)) ;

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

M = M' ;
Ve = Ve' ;
