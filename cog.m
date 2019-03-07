%IMPORT FLIGHT DATA FILE AND FUEL MOMENTS
%Cit_par_mat;
%statespace_EOM;
%FD = open('FLIGHTDATA.m');
FUELfile = csvread('FUELdata.txt') ;

%Static weights in POUNDS
BEM       =  8060 ;    %Provided Basic Empty Mass in POUNDS
FUEL0     =  730 ;    %Total Fuel on T=0
PAYLOAD   =  680 ;    %Total payload in POUNDS
Mass_ramp =  BEM + FUEL0 + PAYLOAD;
ZFM       =  BEM + PAYLOAD;
MAC       =  89.89 ;           %inches

%START CENTER OF GRAVITY
D_pax = [131,131,214,214,251,251,288,288,170];   %inch
D_bag = [74,321,338];                            %inch
Mass_pax = [0,0,0,0,0,0,0,0,0];               %pounds
Mass_bag = [0,0,0];                            %pounds
Moment_pax = D_pax.*Mass_pax;
Moment_bag = D_bag.*Mass_bag;
Moment_pay  = sum(Moment_pax) + sum(Moment_bag);

Weight_nose = 0;
Weight_main = 0;
D_bem = 300.21 - (218.2*Weight_nose / (Weight_nose+Weight_main));
Moment_BEM  = BEM * D_bem;

Moment_FUEL0 = FUEL0 *interp1(FUELfile(:,1),FUELfile(:,2),FUEL0);

Moment_ZFM  = Moment_BEM + Moment_pay;
Moment_RAMP = Moment_ZFM + Moment_FUEL0;
CG_start = Moment_RAMP / Mass_ramp     ;     %defined from nose

%CONTINUOUS MASS AND CG THROUGOUT FLIGHT DATA
%for i=1:length(FD)
%    Mass_t(i)= Mass_ramp - integral(flow,0,t);
%end


%CONVERSION FUNCTIONS
function pounds = kg2p(kg)
    pounds = kg / 0.453592;
end
function kg = p2kg(pounds)
    kg = pounds * 0.453592;
end
function inch = cm2inch(cm)
    inch = cm/2.54;
end
function cm = inch2cm(inch)
    cm = inch*2.54;
end


    