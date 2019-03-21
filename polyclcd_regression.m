% %% import data
% flightdata = open('flightdata.mat')
% Cit_par
% 
% % constants
% g = 9.81 ;          % [m/s^2] gravitational constant
% gamma = 1.4 ;       % spec. heat ratio
% R = 287.05 ;        % [J/kg*K] - spec. gas constant
% S = 30.00 ;	        % [m^2] wing area
% c = 2.0569;         % mean aerodynamic cord [m]
% E_d = 0.686 ;       % [m] characteristic diameter JT15D-4 engine
% bpr = 2.6 ;         % by-pass ratio
% mf = 35.245 ;       % [kg/s] - total mass flow
% mff = 0.177;        % [kg/s] - engine fuel flow
% mffs = 0.048;       % [kg/s] - standard fuel flow per engine
% lambda = -0.0065;   % [K/m] - temperature gradient in ISA
% Ws = 60500 ;        % [N] - standard aircraft weight
% Cmde = -1.1642 ;    % Elevator effectiveness
% CmTc = -0.0064 ;    % dimensionless thrust moment arm
% W_op = 5000*g; % N
% 
% % sea-level ISA parameters
% [T_isa, a_isa, P_isa, rho_isa] = atmoscoesa(0) ;
% 
% % series 1 measurements
% % t = [18*60+3, 20*60+4, 24*60+1, 27*60+26, 29*60+37, 37*60+0]'; % sec
% % hp0 = [6990, 7000, 8010, 8000, 8010, 8980]'; % ft
% % Vc = [241, 219, 192, 161, 137, 119]'; % kts
% % TAT = [1.2, -1.2, -2.2, -6.5, -6.2, -7.0]'; % deg C
% % mfl = [720, 603, 492, 464, 387, 373]'; % lbs/hr
% % mfr = [756, 643, 533, 509, 426, 409]'; % lbs/hr
% % Wfu = [370, 407, 481, 530, 559, 670]'; % lbs
% % alpha = [1.6, 2.4, 3.5, 5.3, 7.7, 10.2]'; % deg
% 
% % reference data
% % t = [19*60+17; 21*60+37; 23*60+46; 26*60+4; 29*60+47; 32*60]; % sec
% % hp0 = [5010; 5020; 5020; 5030; 5020; 5110]; % ft
% % Vc = [249; 221; 192; 163; 130; 118]; % kts
% % TAT = [12.5; 10.5; 8.8; 7.2; 6; 5.2]; % deg C
% % mfl = [798; 673; 561; 463; 443; 474]; % lbs/hr
% % mfr = [813; 682; 579; 484; 467; 499]; % lbs/hr
% % Wfu = [360; 412; 447; 478; 532; 570]; % lbs
% % alpha = [1.7; 2.4; 3.6; 5.4; 8.7; 10.6]; % deg

% convert units
% measured total temp.
TAT = TAT + 273.15 ;    % [K]
% pressure altitude
hp0 = hp0./3.2808399 ;  % [m]
% fuel flow left engine
mfl = mfl.*(0.45359237/3600) ;   % [kg/s]
% fuel flow right engine
mfr = mfr.*(0.45359237/3600) ; % [kg/s]
% fuel used
Wfu = Wfu *0.45359237; % [kg]
% calibrated airspeed 
Vc = Vc.*0.51444 ; %[m/s]

%% Get data points
[Mach, Ve, rho] = speedred(Vc, TAT, hp0, T_isa, P_isa, rho_isa, g, R, gamma, lambda) ;

[T, Ts] = thrustcalc(TAT, hp0, lambda, Mach, mfl, mfr, mffs) ;

[Mass, W, Ve_r, Mass_ramp] = massred(Wfu, Ws, Ve, g) ;

[CL, CD] = liftdrag_coeff(W, T, Ve, S);


%% regression CD = CD0 + k* CL^2
[c2, c1, cinf] = regression(CL, CD);
[CD0_flight, e_flight] = oswald(c2, c1, cinf, A);
CL_alpha = CLalpha(CL, alpha);

%% Plot result
% 
% % plot data points
% plot(CL,CD, 'o')
% %xlim([x(1) - 1, x(m) + 1]) % set x-limits from 0 to 8
% hold on % new plots will join in this window instead of overwriting 
% x = 0 :0.01: CL(6)*1.2; % set a range of x values for plotting
% % calculate range of y values for plotting of regression curves
% yinf = cinf(1) + cinf(2)*x.^2;
% y1 = c1(1) + c1(2)*x.^2;
% y2 = c2(1) + c2(2)*x.^2;
% ygiven = CD0 + x.^2 /(pi*A*e);
% plot(x,y2) % plot regression curves
% plot(x, y1);
% plot(x, yinf);
% plot(x, ygiven)
% % add a legend
% legend({'data points','2-norm optimization', '1-norm optimization', ...
%     'infinity norm optimization', 'given'},'Location','northwest')
% 
% % CL-alpha plot
% figure
% plot(alpha, CL)

%% functions
% P2KG
% KG2P
% liftdrag_coeff
% oswald
% CLalpha
% regression
% speedred
% thrustcalc
% massred

function pounds = KG2P(kg)
    pounds = kg / 0.453592;    
end

function kg = P2KG(pounds)
    kg = pounds * 0.453592;    
end

function [CL, CD] = liftdrag_coeff(W, T, Ve, S)
CL = 2*W./(1.225*Ve.^2*S);
CD = 2*T./(1.225*Ve.^2*S);
end

function [CD0_flight, e_flight] = oswald(c2, c1, cinf, A)
CD0_flight = [c2(1); c1(1); cinf(1)];
e_flight = [c2(2); c1(2); cinf(2)];
e_flight = 1./(pi*A*e_flight);
end

function CL_alpha = CLalpha(CL, alpha)
m = max(size(CL)); %nr of data points
X = [ones(m,1), alpha];
CL_alpha = ((transpose(X)*X)^(-1))*transpose(X) * CL;
CL_alpha = CL_alpha(2);
end


function [c2, c1, cinf] = regression(CL, CD)
m = max(size(CL)); %nr of data points

X = [ones(m,1), CL.^2];

% closed for optimal least squares solution
% copt = ((transpose(X)*X)^(-1))*transpose(X) * CD;

% cvx optimization
cvx_begin quiet
variable c(2); % variable is a vector with 4 elements
minimize norm((X*c-CD),2); %function to minimize. Note that we are not taking the
% norm squared. A norm is always positive so minimizing this is the same
% problem as optimizing the norm squared.
cvx_end
c2 = c;

cvx_begin quiet
variable c(2); % variable is a vector with 4 elements
minimize norm((X*c-CD),1); %function to minimize. Note that we are not taking the
% norm squared. A norm is always positive so minimizing this is the same
% problem as optimizing the norm squared.
cvx_end
c1 = c;

cvx_begin quiet
variable c(2); % variable is a vector with 4 elements
minimize norm((X*c-CD),inf); %function to minimize. Note that we are not taking the
% norm squared. A norm is always positive so minimizing this is the same
% problem as optimizing the norm squared.
cvx_end
cinf = c;
end

function [Mach, Ve, rho] = speedred(Vc, TAT, hp0, T_isa, P_isa, rho_isa, g, R, gamma, lambda)

    for i=1:max(size(Vc))

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
    
    % standard fuel flow list
    mffs = ones(max(size(TAT)), 1) * mffs;

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

    for j=1:max(size(Wfu))

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
    
    Mass_ramp = P2KG(Mass_ramp);
    Mass = Mass' ;
    W = W' ;
    Ve_r = Ve_r' ;
    
end