load flightdata.mat

% ---------------------------------------------------------------------

% measured total temp.
Tm_t = TAT ; % [C]
Tm_t = Tm_t + 273.15 ; % [K]

% pressure altitude
h_p = h; % [ft]
h_p = h_p./3.2808399 ; % [m]

% fuel flow left engine
mfl = FFl ;  % [lbs/h]
mfl = mfl*(0.45359237/3600) ; % [kg/s]

% fuel flow right engine
mfr = FFr ;  % [lbs/h]
mfr = mfr*(0.45359237/3600) ; % [kg/s]

% Mach number
M = flightdata.Dadc1_mach.data ; % [M]
for i=1:6
g_0 = g ; % gravitational constant
gamma = 1.4 ;   % spec. heat ratio
R_s = 287.058 ; % [J/kg*K] - spec. gas constant
bpr = 2.6 ;     % by-pass ratio
m_f1 = FFl(i) ;  % [kg/s] - fuel flow left engine
m_f2 = FFr(i) ;  % [kg/s] - fuel flow right engine
V_c = V(i) ; % [m/s] - calibrated airspeed
T_m = TAT(i) ;     % [K] - measured air temp.

h_p = h(i);   % [m] - ISA pressure altitude
[T_0, a_0, P_0, rho_0] = atmoscoesa(h_p) ; % ISA parameters

% ---------- Reduction of the airspeed -------------

% static pressure
P = P_0*(1 + ((bpr*h_p)/T_0))^(g_0/(bpr*R_s)) ;

% Mach number
M = sqrt((2/(gamma - 1))*((1 + (P_0/P)*((1 + ((gamma - 1)*rho_0*(V_c)^2)/(2*gamma*P_0))^(gamma/(gamma - 1)) - 1))^((gamma - 1)/gamma) - 1)) ;

M(i) = M;
end
% ISA Temperature
T_ISA = 288.15 + ((-0.0065)*h_p) ; % [K]

% Temp. difference
dT = Tm_t - T_ISA ; % [K]

% ---------------------------------------------------------------------

% Merging data set - correct format
data = [h_p M dT mfl mfr] ;
dlmwrite('matlab.dat',data,'delimiter',' ') 
load matlab.dat   % load the file