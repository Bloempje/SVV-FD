% ------------------- Import Flight data  --------------------

FD = open('FTISxprt-20190308_125059.mat') ;     
Time = FD.flightdata.time.data ;    %Time in (mili)seconds

% ----------------------- Constants --------------------------

g = 9.81 ;          % [m/s^2] gravitational constant
gamma = 1.4 ;       % spec. heat ratio
R = 287.05 ;        % [J/kg*K] - spec. gas constant
bpr = 2.6 ;         % by-pass ratio
mf = 35.245 ;       % [kg/s] - total mass flow
mff = 0.177;        % [kg/s] - engine fuel flow
mffs = 0.048 ;      % [kg/s] - standard engine fuel flow per engine
lambda = -0.0065;   % [K/m] - temperature gradient in ISA
Ws = 60500 ;        % [N] - standard aircraft weight
Cmde = -1.1642 ;    % Elevator effectiveness
CmTc = -0.0064 ;    % dimensionless thrust moment arm

% -------------------------------------------------------------------

% datasheet
hp0 = [11920, 12050, 12180, 12270, 11310, 10750, 10280]' ;  % [ft]
Vc = [160, 149, 138, 130, 170, 181, 190]' ;                 % [kts]
TAT = [-11.8, -12.8, -13.8, -14.2, -10, -8.2, -6.2]' ;      % [C]
mfl = [391, 389, 388, 385, 402, 410, 417]' ;                % [lbs/hr]
mfr = [421, 419, 419, 416, 435, 444, 451]' ;                % [lbs/hr]
Wfu = [762, 778, 802, 817, 839, 854, 870]' ;                % [lbs]
de = [0.2, -0.2, -0.7, -1.1, 0.5, 0.8, 1.1]' ;              % [degrees]

% sea-level ISA parameters
[T_isa, a_isa, P_isa, rho_isa] = atmoscoesa(0) ;

% ----------------- Conversion of units ---------------------

% measured total temp.
TAT = TAT + 273.15 ;    % [K]

% pressure altitude
hp0 = hp0./3.2808399 ;  % [m]

% fuel flow left engine
mfl = mfl.*(0.45359237/3600) ;   % [kg/s]

% fuel flow right engine
mfr = mfr.*(0.45359237/3600) ;   % [kg/s]

% ---------- Reduction of the measured airspeed -------------

for i=1:7
    
    V_c = Vc(i)*0.51444 ;   % [m/s] - calibrated airspeed
    T_m = TAT(i) ;          % [K] - measured air temp.
    h_p = hp0(i);           % [m] - ISA pressure altitude

    % static pressure
    P = P_isa*(1 + ((lambda*h_p)/T_isa))^(g/(lambda*R)) ;

    % Mach number
    Mch = sqrt((2/(gamma - 1))*((1 + (P_isa/P)*((1 + ((gamma - 1)*rho_isa*(V_c)^2)/(2*gamma*P_isa))^(gamma/(gamma - 1)) - 1))^((gamma - 1)/gamma) - 1)) ;

    Mj(i) = Mch ;
    
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

M = Mj' ;
Ve = Ve' ;

% --------------------- Thrust calculations --------------------------

% ISA Temperature
T_ISA = 288.15 + (lambda*hp0) ; % [K]

% Temp. difference
dT = TAT - T_ISA ;              % [K]

% Merging data set - correct format
data = [hp0 M dT mfl mfr] ;
data = round(data*1000)/1000 ;
dlmwrite('matlab.dat',data,'delimiter',' ') 
load matlab.dat   % load the file

% Execute thrust calculations
!Thrust.exe &

pause(1) ;

load thrust.dat   % load thrust force per engine

% ---------- Reduction of the non-standard aircraft mass -------------

%START MASS IN POUNDS
BEM  =  9165 ;    %Provided Basic Empty Mass in POUNDS
FUEL0=  3850 ;    %Total Fuel on T=0

Mass_pax  = KG2P([82,92,56,63,75,75,76,77,75]) ; 
Mass_bag  = KG2P([0,0,0]) ;                       
Mass_PAY  = sum(Mass_pax) + sum(Mass_bag) ;   %Total payload in POUNDS

Mass_ramp =  BEM + FUEL0 + Mass_PAY ;
ZFM       =  BEM + Mass_PAY ;

for i=1:7
    
    % aircraft mass
    Mass_t = Mass_ramp - Wfu(i) ; % [lbs]
    
    Mass_t = P2KG(Mass_t) ; % [kg]
    
    Mass(i) = Mass_t ;
    
    % aircraft weight
    Wght = Mass_t*g ;
    
    W(i) = Wght;    % [N]
    
    % reduced equivalent airspeed
    Ver = Ve(i)*sqrt(Ws/Wght) ; %[m/s]
    Ve_r(i) = Ver ;
    
end

Mass = Mass' ;
W = W' ;
Ve_r = Ve_r' ;

% ---------- Reduction of the non-standard engine thrust -------------



% ----------------- Conversion functions --------------------

function pounds = KG2P(kg)
    pounds = kg / 0.453592;
end
function kg = P2KG(pounds)
    kg = pounds * 0.453592;
end


