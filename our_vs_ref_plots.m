% import data
flightdata = open('flightdata.mat')
Cit_par

% constants
g = 9.81 ;          % [m/s^2] gravitational constant
gamma = 1.4 ;       % spec. heat ratio
R = 287.05 ;        % [J/kg*K] - spec. gas constant
S = 30.00 ;	        % [m^2] wing area
c = 2.0569;         % mean aerodynamic cord [m]
E_d = 0.686 ;       % [m] characteristic diameter JT15D-4 engine
bpr = 2.6 ;         % by-pass ratio
mf = 35.245 ;       % [kg/s] - total mass flow
mff = 0.177;        % [kg/s] - engine fuel flow
mffs = 0.048;       % [kg/s] - standard fuel flow per engine
lambda = -0.0065;   % [K/m] - temperature gradient in ISA
Ws = 60500 ;        % [N] - standard aircraft weight
Cmde = -1.1642 ;    % Elevator effectiveness
CmTc = -0.0064 ;    % dimensionless thrust moment arm
W_ramp = 6574.5*g;  % [N]

% sea-level ISA parameters
[T_isa, a_isa, P_isa, rho_isa] = atmoscoesa(0) ;

% series 1 measurements
t = [18*60+3, 20*60+4, 24*60+1, 27*60+26, 29*60+37, 37*60+0]'; % sec
hp0 = [6990, 7000, 8010, 8000, 8010, 8980]'; % ft
Vc = [241, 219, 192, 161, 137, 119]'; % kts
TAT = [1.2, -1.2, -2.2, -6.5, -6.2, -7.0]'; % deg C
mfl = [720, 603, 492, 464, 387, 373]'; % lbs/hr
mfr = [756, 643, 533, 509, 426, 409]'; % lbs/hr
Wfu = [370, 407, 481, 530, 559, 670]'; % lbs
alpha = [1.6, 2.4, 3.5, 5.3, 7.7, 10.2]'; % deg

polyclcd_regression
CL_our = CL;
CD_our = CD;
alpha_our = alpha;
copt_our = c2

% reference data
t = [19*60+17; 21*60+37; 23*60+46; 26*60+4; 29*60+47; 32*60]; % sec
hp0 = [5010; 5020; 5020; 5030; 5020; 5110]; % ft
Vc = [249; 221; 192; 163; 130; 118]; % kts
TAT = [12.5; 10.5; 8.8; 7.2; 6; 5.2]; % deg C
mfl = [798; 673; 561; 463; 443; 474]; % lbs/hr
mfr = [813; 682; 579; 484; 467; 499]; % lbs/hr
Wfu = [360; 412; 447; 478; 532; 570]; % lbs
alpha = [1.7; 2.4; 3.6; 5.4; 8.7; 10.6]; % deg

polyclcd_regression
CL_ref = CL;
CD_ref = CD;
alpha_ref = alpha;
copt_ref = c2
% CL_ref = [0.2148    0.2722    0.3601    0.4990    0.7833    0.9502];
% CD_ref = [0.0242    0.0251    0.0271    0.0297    0.0480    0.0678];
% copt_ref = [0.0202
%     0.0499];
% alpha_ref = [1.7000
%     2.4000
%     3.6000
%     5.4000
%     8.7000
%    10.6000];
% 
% CL_our = [0.2280    0.2779    0.3611    0.5123    0.7064    0.9351];
% CD_our = [0.0237    0.0239    0.0247    0.0352    0.0385    0.0516];
% copt_our = [0.0221
%     0.0342];
% alpha_our = [1.7000
%     2.4000
%     3.5000
%     5.3000
%     7.6000
%    10.3000];

%% Lift-Drag Polar plot
figure
plot(CL_our,CD_our, 'bo')
hold on
plot(CL_ref,CD_ref, 'r*')
%xlim([x(1) - 1, x(m) + 1]) % set x-limits from 0 to 8
hold on % new plots will join in this window instead of overwriting 
x = 0 :0.01: 1.2; % set a range of x values for plotting
% calculate range of y values for plotting of regression curves
y_our = copt_our(1) + copt_our(2)*x.^2;
y_ref = copt_ref(1) + copt_ref(2)*x.^2;
plot(x,y_our) % plot regression curves
plot(x, y_ref, '--');
% title('Lift-Drag Polar')
xlabel('CL')
ylabel('CD')
text(0.05, 0.075 ,'Mach = [0.2-0.4]','FontSize',14)
text(0.05, 0.07 ,'Re = [8 - 16] e6','FontSize',14)
% add a legend
legend({'flight test data points','reference data points', ...
    'flight test regression', 'reference data regression'},...
    'Location','northwest')
%% CL-alpha plot
figure 
plot(alpha_our, CL_our); hold on
plot (alpha_ref, CL_ref, '--')
xlabel('Angle of Attack')
ylabel('CL')
% title('CL-alpha')
text(1.3, 0.8 ,'Mach = [0.2-0.4]','FontSize',14)
text(1.3, 0.75 ,'Re = [8 - 16] e6','FontSize',14)
legend({'flight test', 'reference data'}, 'Location' , 'northwest')