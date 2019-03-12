g_0 = 9.80665 ; % gravitational constant
gamma = 1.4 ;   % spec. heat ratio
R_s = 287.058 ; % [J/kg*K] - spec. gas constant
bpr = 2.6 ;     % by-pass ratio
m_f1 = 0.177 ;  % [kg/s] - fuel flow left engine
m_f2 = 0.177 ;  % [kg/s] - fuel flow right engine
m = 35.245 ;    % [kg/s] - total mass flow
V_c = 160/3.6 ; % [m/s] - calibrated airspeed
T_m = 200 ;     % [K] - measured air temp.

h_p = 2000 ;   % [m] - ISA pressure altitude
[T_0, a_0, P_0, rho_0] = atmoscoesa(h_p) ; % ISA parameters

% ---------- Reduction of the airspeed -------------

% static pressure
P = P_0*(1 + ((bpr*h_p)/T_0))^(g_0/(bpr*R_s)) ;

% Mach number
M = sqrt((2/(gamma - 1))*((1 + (P_0/P)*((1 + ((gamma - 1)*rho_0*(V_c)^2)/(2*gamma*P_0))^(gamma/(gamma - 1)) - 1))^((gamma - 1)/gamma) - 1)) ;

% correct for ram rise
T = T_m/(1+((gamma - 1)/2)*M^2) ;

% density based on ideal gas law
rho_m = P/(R_s*T) ;

% local speed of sound
a = sqrt(gamma*R_s*T) ;

% true airspeed
V_t = M*a ;
% equivalent airspeed
V_e = V_t*sqrt(rho_m/rho_0) ;