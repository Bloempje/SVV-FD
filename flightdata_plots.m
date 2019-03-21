

% select the data from the fligth test validation data 
%phugoid
[time_ph, u_ph,AoA_ph, theta_ph, q_ph] = select_data(time, t1_phugoid, t2_phugoid, true_V, AoA, pitch, q)
%short period
[time_sp, u_sp,AoA_sp, theta_sp, q_sp] = select_data(time, t1_shortperiod, t2_shortperiod, true_V, AoA, pitch, q)


%dutch roll
[time_dr, B_dr,roll_dr, p_dr, r_dr] = select_data(time, t1_dutchroll, t2_dutchroll, sideslip ,roll_angle, p, r)
%spiral
[time_spiral, B_spiral,roll_spiral, p_spiral, r_spiral] = select_data(time, t1_spiral, t2_spiral, sideslip ,roll_angle, p, r)

%aperiodic roll
[time_roll, B_roll,roll_roll, p_roll, r_roll] = select_data(time, t1_aperroll, t2_aperroll, sideslip ,roll_angle, p, r)


%creating time arrays starting at 0
time_ph_plot=time_ph-time_ph(1)
time_sp_plot=time_sp-time_sp(1)
time_dr_plot=time_dr-time_dr(1)
time_spiral_plot=time_spiral-time_spiral(1)
time_roll_plot=time_roll-time_roll(1)


Period_test_dutchroll= Period_test_data(p_roll)
Period_test_phugoid= Period_test_data(u_ph)
Period_test_shortperiod= Period_test_data(AoA_sp)

%plotting the simulation data together with the validation data 
figure
subplot(2,2,1)
plot(time_ph_plot,u_ph,time_ph_plot, phugoid(:,1)+V_ph)
xlabel('Elapsed Time [s]')
ylabel('Velocity [m/s]')

legend({'test', 'simulation'})
subplot(2,2,2)
plot(time_ph_plot,AoA_ph,time_ph_plot,rad2deg(phugoid(:,2))+a0_ph)
xlabel('Elapsed Time [s]')
ylabel('Angle Of Attack [deg]')
axis([time_ph_plot(1) time_ph_plot(end) 0 10])

subplot(2,2,3)
plot(time_ph_plot,theta_ph,time_ph_plot,rad2deg(phugoid(:,3))+th0_ph)
xlabel('Elapsed Time [s]')
ylabel('Pitch Angle [deg]')

subplot(2,2,4)
plot(time_ph_plot,q_ph,time_ph_plot,rad2deg(phugoid(:,4)))
xlabel('Elapsed Time [s]')
ylabel('Pitch rate [deg/s]')




figure
subplot(2,2,1)
plot(time_sp_plot,u_sp,time_sp_plot, shortperiod(:,1)+V_sp)
legend({'test', 'simulation'})
xlabel('Elapsed Time [s]')
ylabel('Velocity [m/s]')
title('velocity')
subplot(2,2,2)
plot(time_sp_plot,AoA_sp,time_sp_plot,rad2deg(shortperiod(:,2))+a0_sp)
xlabel('Elapsed Time [s]')
ylabel('Angle of Attack [deg]')
subplot(2,2,3)
plot(time_sp_plot,theta_sp,time_sp_plot,rad2deg(shortperiod(:,3))+th0_sp)
xlabel('Elapsed Time [s]')
ylabel('Pitch Angle [deg]')
subplot(2,2,4)
plot(time_sp_plot,q_sp,time_sp_plot,rad2deg(shortperiod(:,4)))
xlabel('Elapsed Time [s]')
ylabel('Pitch Rate [deg/s]')

figure
subplot(2,1,2)
plot(time_sp_plot,rad2deg(input_sp))
xlabel('Elapsed Time [s]')
ylabel('Elevator Deflection[deg]')
title(' Short Period - Control input')
subplot(2,1,1)
plot(time_ph_plot,rad2deg(input_ph))
xlabel('Elapsed Time [s]')
ylabel('Elevator Deflection[deg]')
title('Phugoid - Control input')



figure
subplot(3,1,1)
plot(time_dr_plot, roll_dr,time_dr_plot,rad2deg(dutchroll(:,2))+roll_angle0_dr)
xlabel('Elapsed Time [s]')
ylabel('Roll Angle[deg]')
title('Dutch Roll')
legend({'test', 'simulation'})
subplot(3,1,2)
plot(time_dr_plot,p_dr,time_dr_plot,rad2deg(dutchroll(:,3))+p0_dr)
xlabel('Elapsed Time [s]')
ylabel('Roll Rate[deg/s]')
subplot(3,1,3)
plot(time_dr_plot, r_dr,time_dr_plot,rad2deg(dutchroll(:,4))+r0_dr)
xlabel('Elapsed Time [s]')
ylabel('Yaw Rate[deg/s]')

figure

subplot(3,1,1)
plot(time_spiral_plot, roll_spiral,time_spiral_plot,rad2deg( spiral(:,2))+roll_angle0_spiral)
xlabel('Elapsed Time [s]')
ylabel('Roll Angle[deg]')
title('Spiral')
legend({'test', 'simulation'})
subplot(3,1,2)
plot(time_spiral_plot, p_spiral,time_spiral_plot,rad2deg( spiral(:,3))+p0_spiral)
xlabel('Elapsed Time [s]')
ylabel('Roll Rate[deg/s]')
subplot(3,1,3)
plot(time_spiral_plot, r_spiral,time_spiral_plot, rad2deg(spiral(:,4))+r0_spiral)
xlabel('Elapsed Time [s]')
ylabel('Yaw Rate[deg/s]')


figure
subplot(3,1,1)
plot(time_roll_plot, roll_roll,time_roll_plot,rad2deg( aperroll(:,2))+roll_angle0_roll)
xlabel('Elapsed Time [s]')
ylabel('Roll Angle [deg]')
title('Aperiodic Roll')
legend({'test', 'simulation'})
subplot(3,1,2)
plot(time_roll_plot, p_roll, time_roll_plot,rad2deg(aperroll(:,3))+p0_roll)
xlabel('Elapsed Time [s]')
ylabel('Roll Rate[deg/s]')
subplot(3,1,3)
plot(time_roll_plot, r_roll,time_roll_plot,rad2deg( aperroll(:,4))+r0_roll)
xlabel('Elapsed Time [s]')
ylabel('Yaw Rate[deg/s]')

figure
subplot(3,1,1)
plot(time_dr_plot,rad2deg(input1_dr),time_dr_plot,rad2deg(input2_dr))
xlabel('Elapsed Time [s]')
ylabel('Control Surface Deflection [deg] ')
title('Dutch roll - control inputs')
legend({'aileron ', 'rudder '})
subplot(3,1,2)
plot(time_spiral_plot,rad2deg(input1_spiral),time_spiral_plot,rad2deg(input2_spiral))
xlabel('Elapsed Time [s]')
ylabel('Control Surface Deflection [deg] ')
title('Spiral - control inputs')
subplot(3,1,3)
plot(time_roll_plot,rad2deg(input1_roll),time_roll_plot,rad2deg(input2_roll))
xlabel('Elapsed Time [s]')
ylabel('Control Surface Deflection [deg] ')
title('Aperiodic Roll- control inputs')



function  [time_motion, motion_1, motion_2, motion_3, motion_4]= select_data(time_array, timemotion1, timemotion2,data1, data2, data3, data4 )
function  index=select_index(time_motion,time)
index=find(time==time_motion)
end
index1=select_index(timemotion1,time_array)
index2=select_index(timemotion2,time_array)
motion_1= data1(index1:index2)
motion_2= data2(index1:index2)
motion_3= data3(index1:index2)
motion_4= data4(index1:index2)
time_motion=time_array(index1:index2)
end 
function period=Period_test_data(motion)
 [~,peaklocs] = findpeaks(motion);
 period = mean(diff(peaklocs));
end



