%% MASS & Center of Gravity Balance MATLAB FILE B15
tic , Cit_par_mat                                                           %start timer
statespace_EOM, %reduction                                                   %open group files                                                                                   
%% IMPORTS: MATLAB FILES // FLIGHT DATA FILE // FUEL DATA
REFFD = load('ReferenceFD.mat') ;                                           %open reference flightdata
FF = csvread('FUELdata.txt') ;                                              %open fuel moment file
FD = open('FTISxprt-20190308_125059.mat') ;                                 %open flightdata file
FU = FD.flightdata.lh_engine_FU.data + FD.flightdata.rh_engine_FU.data ;    %FUEL USED

Time = FD.flightdata.time.data ;                                            %Time in (deci)sec array
TAS = FD.flightdata.Dadc1_tas.data ;                                        %True Airspeed array
stickforce = FD.flightdata.column_fe.data ;                                 %Stickforce array
AOA = FD.flightdata.vane_AOA.data ;                                         %Angle of Attack array
H_BC = FD.flightdata.Dadc1_bcAlt.data ;                                     %Baro corr. alt array
%% RAMP T=0 START MASS & BALANCE SHEET
BEM        =  9165 ;                                                        %Provided Basic Empty Mass in (POUNDS)
Mass_pax   =  KG2P([82,92,56,63,75,75,76,77,75]);                           %Passenger Mass (INPUT = KG)
Mass_bag   =  KG2P([0,0,0]);                                                %Bags Mass (INPUT = KG)
Mass_pay   =  sum(Mass_pax) + sum(Mass_bag) ;                               %Total payload in POUNDS
Mass_ZF    =  BEM + Mass_pay;                                               %Zero fuel mass (POUNDS)
Mass_fuel0 =  3850 ;                                                        %Total Fuel on T=0 (POUNDS)   
Mass_0     =  BEM + Mass_pay + Mass_fuel0 ;                                 %Ramp mass T=0 (POUNDS)

%RAMP T=0 START CENTER OF GRAVITY
W_nose_jack = 1220*4.44822 ;                                                %Weight on nose jack (N)
W_main_jack = 9945*4.44822 ;                                                %Weight on main jack (N)
D_bem_jack = 315.5 - (221.8*W_nose_jack / (W_nose_jack+W_main_jack));       %ARM BEM (INCH)
D_BEM = 292.18 ;                                                            %ARM BEM provided (INCH)
Moment_BEM  = BEM * D_BEM;                                                  %Moment BEM (INCH-POUNDS)

D_pax = [131,131,214,214,251,251,288,288,170];                              %ARM Passengers (INCH)
D_bag = [74,321,338];                                                       %ARM Bags (INCH)
Moment_pax = D_pax.*(Mass_pax) ;                                            %Moment Passengers (INCH-POUNDS)
Moment_bag = D_bag.*(Mass_bag);                                             %Moment Bags (INCH-POUNDS)
Moment_pay  = sum(Moment_pax) + sum(Moment_bag);                            %Moment Payload (INCH-POUNDS)

Moment_fuel0 = 100*interp1(FF(:,1),FF(:,2),Mass_fuel0);                     %Moment Fuel (INCH-POUNDS)

Moment_0 = Moment_BEM + Moment_pay + Moment_fuel0;                          %Ramp T=0 Moment (INCH-POUNDS)
CG_0 = Moment_0/Mass_0 ;                                                    %Ramp T=0 Xcog (INCH)

%% CONTINUOUS MASS AND CENTER OF GRAVITY                                    %Function of time
Mass_t=zeros(size(Time)); Mass_fuel_t=zeros(size(Time));                    %Preallocating arrays
Moment_fuel_t=zeros(size(Time)); Moment_Shift=zeros(size(Time));            %Preallocating arrays
Moment_t=zeros(size(Time));CG_t=zeros(size(Time));                          %Preallocating arrays

for i=1:length(Time)
    Mass_t(i)= Mass_0 - FU(i) ;                                             %Continuous Mass
    Mass_fuel_t(i) = Mass_fuel0 - FU(i) ;                                   %Continuous Fuel Mass
    Moment_fuel_t(i)= 100*interp1(FF(:,1),FF(:,2),Mass_fuel_t(i));          %Continuous Fuel Moment
    if (i > ((52.5*60*10)-89)) && (i < ((54*60*10)-89))                     %IF time in range Pax shift
        Moment_Shift(i) = Mass_pax(8)*D_pax(1) - Mass_pax(8)*D_pax(8) ;     %Passenger Shift
    else    
        Moment_Shift(i) = 0 ;
    end
    Moment_t(i) = Moment_BEM + Moment_fuel_t(i) + Moment_pay + Moment_Shift(i); %Continuous Moment
    CG_t(i) = Moment_t(i)/Mass_t(i) ;                                           %Continuous Xcg 
end   

%% PLOTS
Line_BEM=zeros(size(Time));Line_PAY=zeros(size(Time));Line_0=zeros(size(Time));
Line_CG1=zeros(size(Time));Line_CG2=zeros(size(Time)) ;                     %preallocating array lines
for i=1:length(Time)
    Line_BEM(i) = BEM ;
    Line_PAY(i) = BEM + Mass_pay;
    Line_0(i) = BEM + Mass_pay + Mass_fuel0;
    Line_CG1(i) = 276.10 ;
    Line_CG2(i) = 285.8 ;
end

figure(21)
plot(Time,Mass_t,Time,FU,Time,Line_BEM,Time,Line_PAY,Time,Line_0);
title('Total Mass & Fuel mass consumption')
legend('Gross Mass','Fuel consumed','BEM','Payload') ;
xlabel('Time(sec)') ;
ylabel('Mass(lbs)') ;
grid on

figure(22)
plot(Time,stickforce,Time,AOA,Time,H_BC/100) ;
title('Stickforce & AOA & ALT - plot') ;
xlabel('Time(sec)') ;
ylabel('variable(?)') ;
grid on

figure(23)
plot(Time,CG_t,Time,Line_CG1,Time,Line_CG2) ;
title('CG plot') ;
xlabel('Time(sec)') ;
ylabel('CG(in)') ;
grid on


%% SI UNITS & EXPORTS
SI_Mass_t = P2KG(Mass_t) ;                                                  %Continuous mass in Kg
SI_CG_t = INCH2CM(CG_t)/100 ;                                               %Continuous cog in M
SI_Moment_t = IP2NM(Moment_t) ;                                             %continuous Moment in Nm

LEMAC = 261.56 ; MAC = 80.89 ;                                              %Leading Edge Mean Aerodynamic Chord                                     
SI_MAC = 2.0569 ;                                                           %MAC in M
SI_CG_MAC = SI_CG_t - LEMAC ;                                               %continuous cog from Lemac in M
SI_CG_PERCHORD = (SI_CG_MAC/c)*100 ;                                        %continuous cog % of chord

CG_t(1) ;
CG_0 ;
CG_t(60121) ;
(Moment_BEM + Moment_pay) / Mass_ZF ;

Mass_t(T2P(50.9))
CG_t(T2P(50.9))

Mass_t(T2P(49.2))
CG_t(T2P(49.2))

toc                                                                         %Stop Timer
%% CONVERSON FUNCTIONS

function datapoint = T2P(minutes)                                           %datapoint -> Time(min)
    datapoint = (minutes*60*10)-89 ; 
end
function minutes = P2T(datapoint)                                           %Time(min) -> datapoint
    minutes = (datapoint + 89)/600 ; 
end
function pounds = KG2P(kg)                                                  %Kg    -> Lbs      
    pounds = kg / 0.453592;     
end
function kg = P2KG(pounds)
    kg = pounds * 0.453592;
end                                               %Lbs   -> Kg
function inch = CM2INCH(cm)
    inch = cm/2.54;
end                                              %Cm    -> Inch
function cm = INCH2CM(inch)
    cm = inch*2.54;
end                                              %Inch  -> Cm
function NewtonMeter = IP2NM(InchPound)
     NewtonMeter = InchPound*0.112984829 ;
end                                  %In-lb -> Nm



    