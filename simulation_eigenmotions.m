%symm input and time array selection
[input_ph, time_ph] = simulation_inputs_sym(time, t1_phugoid, t2_phugoid, deg2rad(elevator))
[input_sp, time_sp] = simulation_inputs_sym(time, t1_shortperiod, t2_shortperiod, deg2rad(elevator))
%symmetrical motion simulation
phugoid= simulation_sym(input_ph, time_ph, system_phugoid, V_ph, a0_ph, th0_ph,q0_ph)
shortperiod=simulation_sym(input_sp, time_sp, system_shortperiod, V_sp,a0_sp,th0_sp ,q0_sp)



%asymmymetric input and time array selection
[input1_dr,input2_dr, timescale_dr]= simulation_inputs_asym(time, t1_dutchroll, t2_dutchroll, deg2rad(aileron), deg2rad(rudder))
[input1_spiral,input2_spiral, timescale_spiral]= simulation_inputs_asym(time, t1_spiral, t2_spiral, deg2rad(aileron), deg2rad(rudder))
[input1_roll,input2_roll, timescale_roll]= simulation_inputs_asym(time, t1_aperroll, t2_aperroll, deg2rad(aileron), deg2rad(rudder))

% assymmetric motion simulation
dutchroll=simulation_asym(input1_dr,input2_dr, timescale_dr, system_dutchroll, 0, roll_angle0_dr,p0_dr,r0_dr)
spiral=simulation_asym(input1_spiral,zeros(length(input1_spiral),1), timescale_spiral, system_spiral, 0, roll_angle0_spiral,p0_spiral,r0_spiral)
aperroll=simulation_asym(input1_roll,input2_roll, timescale_roll, system_roll, 0, roll_angle0_roll,p0_roll,r0_roll)


%functions needed 
function [input1, timescale]= simulation_inputs_sym(time_array, timemotion1, timemotion2, input1_array)
function  index=select_index(time_motion,time)
index=find(time==time_motion)
end
index1=find(time_array==timemotion1)
index2=find(time_array==timemotion2)

timescale=time_array(index1:index2)
input1=input1_array(index1:index2)
input1=input1-input1(1)
%input2=input2_array(index1:index2)
end 

function [input1,input2, timescale]= simulation_inputs_asym(time_array, timemotion1, timemotion2, input1_array,input2_array)
function  index=select_index(time_motion,time)
index=find(time==time_motion)
end
index1=find(time_array==timemotion1)
index2=find(time_array==timemotion2)

timescale=time_array(index1:index2)
input1=input1_array(index1:index2)
input2=input2_array(index1:index2)
input1=input1-input1(1)
input2=input2-input2(1)
end 

function [motion]= simulation_sym( input,timescale, system,  V, a0, th0,q0)
x0=zeros(4,1)
motion=lsim(system, input, timescale,x0)

end

function [motion]= simulation_asym(input1,input2, timescale, system, B0, roll_angle0,p0,r0)
input=[input1 input2]
x0=zeros(4,1)
motion=lsim(system, input, timescale,x0);

end