%timing of all the eigenmotions

t1_phugoid=55*60+5;
t2_phugoid=58*60+0;

t1_shortperiod=59*60+12;
t2_shortperiod=59*60+20;

t1_aperroll=1*3600+2*60+20;
t2_aperroll=1*3600+3*60+15;

t1_dutchroll=1*3600+0*60+15;
t2_dutchroll=1*3600+1*60;

t1_spiral=1*3600+6*60;
t2_spiral=1*3600+8*60;

%opening all the required sets for simulation 

time=flightdata.time.data; %[s]
pressure_alt=flightdata.Dadc1_alt.data *0.3048 %[m]

calibrated_V=flightdata.Dadc1_cas.data*knots_to_m; %[m/s]
true_V=flightdata.Dadc1_tas.data*knots_to_m; %[m/s}

elevator=flightdata.delta_e.data ; %[deg]
stick_force=flightdata.column_fe.data ; %[N]
aileron=-flightdata.delta_a.data; %[deg]
rudder=-flightdata.delta_r.data; %[deg]



%symmetrical
pitch=flightdata.Ahrs1_Pitch.data;   %[deg]
AoA=flightdata.vane_AOA.data;    %[deg]
q=flightdata.Ahrs1_bPitchRate.data ; %[deg/s]


%asymmetrical
sideslip=flightdata.vane_AOA.data;    %[deg]
roll_angle=flightdata.Ahrs1_Roll.data ;   %[deg]
p=flightdata.Ahrs1_bRollRate.data ;   %[deg/s]
r=flightdata.Ahrs1_bYawRate.data ;     %[deg/s]



plot(time, aileron, time, rudder)




