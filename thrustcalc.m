load FTISxprt-20190308_125059.mat

% ---------------------------------------------------------------------

% measured total temp.
Tm_t = flightdata.Dadc1_tat.data ; % [C]
Tm_t = Tm_t + 273.15 ; % [K]

% pressure altitude
h_p = flightdata.Dadc1_alt.data ; % [ft]
h_p = h_p./3.2808399 ; % [m]

% fuel flow left engine
mfl = flightdata.lh_engine_FMF.data ;  % [lbs/h]
mfl = mfl*(0.45359237/3600) ; % [kg/s]

% fuel flow right engine
mfr = flightdata.rh_engine_FMF.data ;  % [lbs/h]
mfr = mfr*(0.45359237/3600) ; % [kg/s]

% Mach number
M = flightdata.Dadc1_mach.data ; % [M]

% ISA Temperature
T_ISA = 288.15 + ((-0.0065)*h_p) ; % [K]

% Temp. difference
dT = Tm_t - T_ISA ; % [K]

% ---------------------------------------------------------------------

% Merging data set - correct format
data = [h_p M dT mfl mfr] ;
data = round(data*1000)/1000 ;
dlmwrite('matlab.dat',data,'delimiter',' ') 
load matlab.dat   % load the file
