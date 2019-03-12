%IMPORT FLIGHT DATA FILE 
FD = open('FTISxprt-20190308_125059.mat') ;     
Time = FD.flightdata.time.data ;                                            %Time in (mili)seconds
FU = FD.flightdata.lh_engine_FU.data + FD.flightdata.rh_engine_FU.data ;    %FUEL USED
TAS = FD.flightdata.Dadc1_tas.data ;                                        %True Airspeed

%IMPORT FUEL DATA AND OTHER MATLAB
%Cit_par_mat , statespace_EOM;
FUELfile = csvread('FUELdata.txt') ;

%START MASS IN POUNDS
BEM  =  9165 ;    %Provided Basic Empty Mass in POUNDS
FUEL0=  3850 ;    %Total Fuel on T=0
Mass_pax  = KG2P([82,92,56,63,75,75,76,77,75]); 
Mass_bag  = KG2P([0,0,0]);                       
Mass_PAY  = sum(Mass_pax) + sum(Mass_bag) ;         %Total payload in POUNDS

Mass_ramp =  BEM + FUEL0 + Mass_PAY ;
ZFM       =  BEM + Mass_PAY;
MAC       =  89.89 ;  
Weight_nose_jack = 1220*9.81;
Weight_main_jack = 9945*9.81;

%START CENTER OF GRAVITY
D_bem_jack = 315.5 - (221.8*Weight_nose_jack / (Weight_nose_jack+Weight_main_jack));
D_BEM = 292.18 ;
D_pax = [131,131,214,214,251,251,288,288,170];   %inch
D_bag = [74,321,338];                            %inch

Moment_BEM  = BEM * D_BEM;
Moment_pax = D_pax.*(Mass_pax);
Moment_bag = D_bag.*(Mass_bag);
Moment_PAY  = sum(Moment_pax) + sum(Moment_bag);
Moment_FUEL0 = 100*interp1(FUELfile(:,1),FUELfile(:,2),FUEL0);

Moment_ZFM  = Moment_BEM + Moment_PAY;
Moment_RAMP = Moment_ZFM + Moment_FUEL0;
CG_start = Moment_RAMP/Mass_ramp ;          %defined from nose

%CONTINUOUS MASS AND
%CONTINUOUS CENTER OF GRAVITY
for i=1:length(FU)
    Mass_t(i)= Mass_ramp - FU(i);
end
for i=1:length(Time)
    BEMline(i) = BEM ;
    BEMPAYline(i) = BEM + Mass_PAY;
    BEMPAYFUELline(i) = BEM + Mass_PAY + FUEL0;
end
plot(Time,Mass_t,Time,FU,Time,BEMplot,Time,BEMPAYplot,Time,BEMPAYFUELline);
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


    