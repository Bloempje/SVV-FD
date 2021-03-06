%% MASS & Center of Gravity Balance MATLAB FILE B15
%% IMPORTS: MATLAB FILES // FLIGHT DATA FILE // FUEL DATA
tic , %Cit_par_mat , statespace_EOM;                                        %open matlab group file                                                                                    
FF = csvread('FUELdata.txt') ;                                              %open fuel moment file
FD = open('FTISxprt-20190308_125059.mat') ;                                 %open flightdata file
Time = FD.flightdata.time.data ;                                            %Time in (deci)seconds
FU = FD.flightdata.lh_engine_FU.data + FD.flightdata.rh_engine_FU.data ;    %FUEL USED

TAS = FD.flightdata.Dadc1_tas.data ;                                        %True Airspeed
stickforce = FD.flightdata.column_fe.data ;                                 %Stickforce
AOA = FD.flightdata.vane_AOA.data ;                                         %Angle of Attack
H_BC = FD.flightdata.Dadc1_bcAlt.data ;                                     %Baro corr. alt
%% START MASS & CG 
BEM  =  9165 ;                                                              %Provided Basic Empty Mass in (POUNDS)
Mass_fuel0=  3850 ;                                                         %Total Fuel on T=0 (POUNDS)   
Mass_pax  = KG2P([82,92,56,63,75,75,76,77,75]);                             %Passenger Mass (INPUT = KG)
Mass_bag  = KG2P([0,0,0]);                                                  %Bags Mass (INPUT = KG)
Mass_pay  = sum(Mass_pax) + sum(Mass_bag) ;                                 %Total payload in POUNDS
Mass_ZF   =  BEM + Mass_pay;                                                %Zero fuel mass (POUNDS)
Mass_0 =  BEM + Mass_pay + Mass_fuel0 ;                                     %Ramp mass T=0 (POUNDS)

%START CENTER OF GRAVITY
W_nose_jack = 1220*4.44822 ;                                                %Weight on nose jack (N)
W_main_jack = 9945*4.44822 ;                                                %Weight on main jack (N)
D_bem_jack = 315.5 - (221.8*W_nose_jack / (W_nose_jack+W_main_jack));       %ARM BEM (INCH)
D_BEM = 292.18 ;                                                            %ARM BEM provided (INCH)
Moment_BEM  = BEM * D_BEM;                                                  %Moment BEM (INCH-POUNDS)

D_pax = [131,131,214,214,251,251,288,288,170];                              %ARM Passengers (INCH)
D_bag = [74,321,338];                                                       %ARM Bags (INCH)
Moment_pax = D_pax.*(Mass_pax);                                             %Moment Passengers (INCH-POUNDS)
Moment_bag = D_bag.*(Mass_bag);                                             %Moment Bags (INCH-POUNDS)
Moment_pay  = sum(Moment_pax) + sum(Moment_bag);                            %Moment Payload (INCH-POUNDS)

Moment_fuel0 = 100*interp1(FF(:,1),FF(:,2),Mass_fuel0);                     %Moment Fuel (INCH-POUNDS)

Moment_0 = Moment_BEM + Moment_pay + Moment_fuel0;                          %Ramp T=0 Moment (INCH-POUNDS)
CG_0 = Moment_0/Mass_0 ;                                                    %Ramp T=0 Xcog (INCH)
%% CONTINUOUS MASS AND  CENTER OF GRAVITY
for i=1:length(Time)
    Mass_t(i)= Mass_0 - FU(i) ;                                             %Continuous Mass
    Mass_fuel_t(i) = Mass_fuel0 - FU(i) ;                                   %Continuous Fuel Mass
    Moment_fuel_t(i)= 100*interp1(FF(:,1),FF(:,2),Mass_fuel_t(i));          %Continuous Fuel Moment
    if (i > ((52.5*60*10)-89)) && (i < ((54.5*60*10)-89))                   %IF time in range Pax shift
        Moment_S(i) = Mass_pax(8)*D_pax(1) - Mass_pax(8)*D_pax(8) ;         %Passenger Shift
    else    
        Moment_S(i) = 0 ;
    end
    Moment_t(i) = Moment_BEM + Moment_fuel_t(i) + Moment_pay + Moment_S(i); %Continuous Moment
    CG_t(i) = Moment_t(i)/Mass_t(i) ;                                       %Continuous Xcog 
end   

%% PLOTS
for i=1:length(Time)
    Line_BEM(i) = BEM ;
    Line_PAY(i) = BEM + Mass_pay;
    Line_0(i) = BEM + Mass_pay + Mass_fuel0;
    Line_CG1(i) = 276.10 ;
    Line_CG2(i) = 285.8 ;
end
%figure
plot(Time,Mass_t,Time,FU,Time,Line_BEM,Time,Line_PAY,Time,Line_0);
title('Total weight & Fuel consumed - plot')
%figure
plot(Time,stickforce,Time,AOA,Time,H_BC/100) ;
title('Stickforce & AOA & ALT - plot') ;
%figure
plot(Time,CG_t,Time,Line_CG1,Time,Line_CG2) ;
title('CG plot') ;

CG_plot = CG_t(1:800:end) ;
Mass_plot = Mass_t(1:800:end)

figure
scatter(CG_plot,Mass_plot,20,[0,0,1])
title('Center of Gravity througout flight')
xlabel('Xcg Fuselage station-inches  (Inch)')
ylabel('Total Mass  (Pounds)')
grid on
grid minor

%% SI UNITS & EXPORTS
SI_Mass_t = P2KG(Mass_t) ; 
SI_CG_t = INCH2CM(CG_t)/100 ; 
SI_Moment_t = IP2NM(Moment_t) ;

LEMAC = INCH2CM(261.56)/100  ;                                              %Leading Edge MAC 
MAC = INCH2CM(80.89)/100     ;                                              %Mean Aerodynamic Chord
c = 2.0569 ;
SI_CG_MAC = SI_CG_t - LEMAC ;
SI_CG_PERCHORD = (SI_CG_MAC/c)*100 ;

SI_CG_PERCHORD(1)
SI_CG_PERCHORD(60121)
SI_CG_PERCHORD((53*60*10)-89)

toc                                                                         %Stop Timer
%% Verification Functions
function som = sommeren(input)
    som = 0 ;
    for i=1:length(input)
        som = som + input(i) ; 
    end
end

function outliers = Listverification(list,ref) 
    outliers = [] ;
    OL = isoutlier(list) ;
    for i=1:length(list) 
        if OL(i) == 1 
        append(outliers,i,list(i)) ; 
        elseif ref(i) ~= 0 
            if sign(list(i)) ~= sign(ref(i)) 
                append(outliers,i,list(i),ref(i)) ;
            elseif floor(log10(list(i))) ~= floor(log10(ref(i)))
                append(outliers,i,list(i),ref(i)) ;
            end 
        end
    end 
end

function change = Crossverification(FUNCTION1,IN1,IN2)
    %Should run MODULE1 with input argument IN1 and IN2
    %outputs of this MODULE1 are the inputs the next MODULE2
    %find change in output of MODULE2 due to IN1 -> IN2
    OUT1 = FUNCTION1(IN1) ; 
    OUT2 = FUNCTION1(IN2) ; 
    disp OUT1
    disp OUT2
end 
%% CONVERSON FUNCTIONS
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
end           