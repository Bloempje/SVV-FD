%IMPORT FLIGHT DATA FILE 
FD = open('FTISxprt-20190308_125059.mat') ;
Time = FD.flightdata.time.data ;
FU = FD.flightdata.lh_engine_FU.data + FD.flightdata.rh_engine_FU.data ;
TAS = FD.flightdata.Dadc1_tas.data ;   

%IMPORT FUEL DATA AND OTHER MATLAB
%Cit_par_mat , statespace_EOM;
FUELfile = csvread('FUELdata.txt') ;

%START MASS IN POUNDS
BEM  =  8060 ;    %Provided Basic Empty Mass in POUNDS
FUEL0=  3850 ;    %Total Fuel on T=0
Mass_pax  = [82,92,75,56,63,75,75,76,77];  %KG!!!
Mass_bag  = [0,0,0];                       %KG!!!
Mass_PAY  =  sum(KG2P(Mass_pax)) ;    %Total payload in POUNDS

Mass_ramp =  BEM + FUEL0 + Mass_PAY;
ZFM       =  BEM + Mass_PAY;
MAC       =  89.89 ;           %inches
Weight_nose = 3000;
Weight_main = 8000;

%START CENTER OF GRAVITY
D_bem = 300.21 - (218.2*Weight_nose / (Weight_nose+Weight_main));
D_pax = [131,131,214,214,251,251,288,288,170];   %inch
D_bag = [74,321,338];                            %inch

Moment_BEM  = BEM * D_bem;
Moment_pax = D_pax.*KG2P(Mass_pax);
Moment_bag = D_bag.*KG2P(Mass_bag);
Moment_PAY  = sum(Moment_pax) + sum(Moment_bag);
Moment_FUEL0 = FUEL0 *interp1(FUELfile(:,1),FUELfile(:,2),FUEL0);

Moment_ZFM  = Moment_BEM + Moment_PAY;
Moment_RAMP = Moment_ZFM + Moment_FUEL0;
CG_start = Moment_RAMP / Mass_ramp     ;     %defined from nose


%CONTINUOUS MASS AND
%CONTINUOUS CENTER OF GRAVITY
for i=1:length(FU)
    Mass_t(i)= Mass_ramp - FU(i);
end

plot(Time,Mass_t,Time,FU);
FDvariables = (fields(FD.flightdata));
    

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


    