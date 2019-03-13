%% MASS & Balance sheet  +  Xcg & Mass shift through flight  
%IMPORTS: MATLAB FILES // FLIGHT DATA FILE // FUEL DATA
Cit_par_mat , statespace_EOM;                                               %open matlab group file
FUELfile = csvread('FUELdata.txt') ;                                        %open fuel moment file
FD = open('FTISxprt-20190308_125059.mat') ;                                 %open flightdata file
Time = FD.flightdata.time.data ;                                            %Time in (mili)seconds
FU = FD.flightdata.lh_engine_FU.data + FD.flightdata.rh_engine_FU.data ;    %FUEL USED
TAS = FD.flightdata.Dadc1_tas.data ;                                        %True Airspeed

%% START MASS & Xcg 
BEM  =  9165 ;    %Provided Basic Empty Mass in POUNDS
Mass_fuel0=  3850 ;    %Total Fuel on T=0
Mass_pax  = KG2P([82,92,56,63,75,75,76,77,75]); 
Mass_bag  = KG2P([0,0,0]);                       
Mass_pay  = sum(Mass_pax) + sum(Mass_bag) ;         %Total payload in POUNDS
Mass_ZF   =  BEM + Mass_pay;
Mass_ramp =  BEM + Mass_pay + Mass_fuel0 ;

%START CENTER OF GRAVITY
Weight_nose_jack = 1220*9.81;
Weight_main_jack = 9945*9.81;
D_bem_jack = 315.5 - (221.8*Weight_nose_jack / (Weight_nose_jack+Weight_main_jack));
D_BEM = 292.18 ;
Moment_BEM  = BEM * D_BEM;

D_pax = [131,131,214,214,251,251,288,288,170];   %inch
D_bag = [74,321,338];                            %inch
Moment_pax = D_pax.*(Mass_pax);
Moment_bag = D_bag.*(Mass_bag);
Moment_pay  = sum(Moment_pax) + sum(Moment_bag);

Moment_fuel0 = 100*interp1(FUELfile(:,1),FUELfile(:,2),Mass_fuel0);

Moment_RAMP = Moment_BEM + Moment_pay + Moment_fuel0;
CG_start = Moment_RAMP/Mass_ramp ;          %defined from nose
%% CONTINUOUS MASS AND  CENTER OF GRAVITY
for i=1:length(FU)
    Mass_t(i)= Mass_ramp - FU(i);
end
for i=1:length(Time)
    Moment_pax_t()
end   
    

%% PLOTS & EXPORTS
for i=1:length(Time)
    BEMline(i) = BEM ;
    BEMpayline(i) = BEM + Mass_pay;
    BEMpayFUELline(i) = BEM + Mass_pay + Mass_fuel0;
end
plot(Time,Mass_t,Time,FU,Time,BEMline,Time,BEMpayline,Time,BEMpayFUELline);
title('Total weight & Fuel consumed - plot')
xlabel('Time (sec)')
ylabel('Weight (lbs)')

FDvariables = (fields(FD.flightdata));    
FU_total = P2KG(FU(length(FU)));



%CONVERSON FUNCTIONS
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


    