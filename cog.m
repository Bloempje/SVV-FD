%% MASS & Center of Gravity Balance MATLAB FILE B15
%% IMPORTS: MATLAB FILES // FLIGHT DATA FILE // FUEL DATA
Cit_par_mat , statespace_EOM;                                               %open matlab group file
FUELfile = csvread('FUELdata.txt') ;                                        %open fuel moment file
FD = open('FTISxprt-20190308_125059.mat') ;                                 %open flightdata file
Time = FD.flightdata.time.data ;                                            %Time in (deci)seconds
FU = FD.flightdata.lh_engine_FU.data + FD.flightdata.rh_engine_FU.data ;    %FUEL USED

TAS = FD.flightdata.Dadc1_tas.data ;                                        %True Airspeed
stickforce = FD.flightdata.column_fe.data ;                                 %Stickforce
AOA = FD.flightdata.vane_AOA.data ;                                         %Angle of Attack
BCalt = FD.flightdata.Dadc1_bcAlt.data ;                                    %Baro corr. alt
%% START MASS & Xcg 
BEM  =  9165 ;                                                              %Provided Basic Empty Mass in (POUNDS)
Mass_fuel0=  3850 ;                                                         %Total Fuel on T=0 (POUNDS)   
Mass_pax  = KG2P([82,92,56,63,75,75,76,77,75]);                             %Passenger Mass (INPUT = KG)
Mass_bag  = KG2P([0,0,0]);                                                  %Bags Mass (INPUT = KG)
Mass_pay  = sum(Mass_pax) + sum(Mass_bag) ;                                 %Total payload in POUNDS
Mass_ZF   =  BEM + Mass_pay;                                                %Zero fuel mass (POUNDS)
Mass_ramp =  BEM + Mass_pay + Mass_fuel0 ;                                  %Ramp mass T=0 (POUNDS)

%START CENTER OF GRAVITY
W_nose_jack = 1220*4,44822 ;                                                %Weight on nose jack (N)
W_main_jack = 9945*4,44822 ;                                                %Weight on main jack (N)
D_bem_jack = 315.5 - (221.8*W_nose_jack / (W_nose_jack+W_main_jack));       %ARM BEM (INCH)
D_BEM = 292.18 ;                                                            %ARM BEM provided (INCH)
Moment_BEM  = BEM * D_BEM;                                                  %Moment BEM (INCH-POUNDS)

D_pax = [131,131,214,214,251,251,288,288,170];                              %ARM Passengers (INCH)
D_bag = [74,321,338];                                                       %ARM Bags (INCH)
Moment_pax = D_pax.*(Mass_pax);                                             %Moment Passengers (INCH-POUNDS)
Moment_bag = D_bag.*(Mass_bag);                                             %Moment Bags (INCH-POUNDS)
Moment_pay  = sum(Moment_pax) + sum(Moment_bag);                            %Moment Payload (INCH-POUNDS)

Moment_fuel0 = 100*interp1(FUELfile(:,1),FUELfile(:,2),Mass_fuel0);         %Moment Fuel (INCH-POUNDS)

Moment_RAMP = Moment_BEM + Moment_pay + Moment_fuel0;                       %Moment Ramp T=0 (INCH-POUNDS)
CG_start = Moment_RAMP/Mass_ramp ;                                          %Xcog (INCH)

%% CONTINUOUS MASS AND  CENTER OF GRAVITY
for i=1:length(Time)
    Mass_t(i)= Mass_ramp - FU(i) ;                                          %
    Mass_fuel_t(i) = Mass_fuel0 - FU(i) ;                                   %
    Moment_fuel_t(i)= 100*interp1(FUELfile(:,1),FUELfile(:,2),Mass_fuel_t(i))  ; %
    if (i > ((52.5*60*10)-89)) && (i < ((54.5*60*10)-89))                   %
        Moment_paxshift(i) = Mass_pax(8)*D_pax(1) - Mass_pax(8)*D_pax(8) ;  %
    else    
        Moment_paxshift(i) = 0 ;
    end
    Moment_t(i) = Moment_BEM + Moment_fuel_t(i) + Moment_pay + Moment_paxshift(i) ;
    CG_t(i) = Moment_t(i)/Mass_t(i) ;
end   

%% PLOTS
for i=1:length(Time)
    BEMline(i) = BEM ;
    BEMpayline(i) = BEM + Mass_pay;
    BEMpayFUELline(i) = BEM + Mass_pay + Mass_fuel0;
end
%plot(Time,Mass_t,Time,FU,Time,BEMline,Time,BEMpayline,Time,BEMpayFUELline);
%title('Total weight & Fuel consumed - plot')

plot(Time,stickforce,Time,AOA,Time,BCalt/1000) ;
title('Stickforce & AOA & ALT - plot') ;

%plot(Time,CG_t) ;
title('Xcg plot') ;
%% SI UNITS EXPORTS
Mass_kg = P2KG(Mass_t) ; 
XCG_cm = INCH2CM(CG_t) ; 


%% CONVERSON FUNCTIONS
function pounds = KG2P(kg)
    pounds = kg / 0.453592;
end
function kg = P2KG(pounds)
    kg = pounds * 0.453592;
end
function inch = CM2INCH(cm)
    inch = cm/2.54;
end
function cm = INCH2CM(inch)
    cm = inch*2.54;
end


    