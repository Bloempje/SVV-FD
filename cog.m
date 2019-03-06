%IMPORT FLIGHT DATA FILE AND FUEL MOMENTS
FD = open('FLIGHTDATA.m')   ;
Timet = FD(:,1)      ;
Fuelt = FD(:,2)      ;
FUELfile = csvread('FUELdata.txt')
FUELround = interp1(FUELfile(:,1),FUELfile(:,2),FUEL0);

%Static weights in POUNDS
BEM       =  8060 ;    %Provided Basic Empty Mass in POUNDS
FUEL0     =  1150 ;    %Total Fuel on T=0
PAYLOAD   =  2350 ;    %Total payload in POUNDS
Mass_ramp = BEM + FUEL + PAYLOAD;
ZFM       = BEM + PAYLOAD
MAC = 89.89            %inches

%START CENTER OF GRAVITY (ACC.  276.1 < INCH < 285.8 FROM NOSE)
D_pax = [131,131,214,214,251,251,288,288,170]   %inch
D_bag = [74,321,338]                            %inch
Mass_pax = [0,0,0,0]                            %pounds
Mass_bag = [0,0,0,0]                            %pounds

Weight_nose = 0
Weight_main = 0
D_bem = 300.21 - (218.2*Weight_nose / (Weight_nose+Weight_main))

for i=1:length(D_pax)                    %INCH POUNDS
    Moment_pax(i) = D_pax(i)*Mass_pax(i)
    Moment_bag(i) = D_bag(i)*Mass_bag(i)
end
    Moment_BEM  = BEM * D_bem
    Moment_ZFM  = Moment_BEM + sum(Moment_pax) + sum(Moment_bag)
    Moment_RAMP = Moment_ZFM + Moment_FUEL0
CG_start = Moment_RAMP / Mass_ramp          %defined from nose

%CONTINUOUS MASS AND CG THROUGOUT FLIGHT DATA
for i=1:length(FD)
    Mass_t(i)= Mass_ramp - integral(flow,0,t);
end


%CONVERSION FUNCTIONS
function pounds = kg2p(kg)
    pounds = kg / 0.453592
end
function kg = p2kg(pounds)
    kg = pounds * 0.453592
end
function inch = cm2inch(cm)
    inch = cm/2.54
end
function cm = inch2cm(inch)
    cm = inch*2.54
end


    