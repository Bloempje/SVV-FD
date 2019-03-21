% INITIAL VALUE RESPONES OF THE NUMERICAL MODEL

x0_symm=[0,0.1,0,0]
dist_sym=initial(system_phugoid,x0_symm)
figure
subplot(4,1,1)
plot(dist_sym(:,1))
subplot(4,1,2)
plot(dist_sym(:,2))
subplot(4,1,3)
plot(dist_sym(:,3))
subplot(4,1,4)
plot(dist_sym(:,4))

x0_asym=[0;0 ; 3 ; 0]
dist_asym=initial(system_dutchroll,x0_asym )

figure
subplot(4,1,1)
plot(dist_asym(:,1))
axis([0 200 0 100])
subplot(4,1,2)
plot(dist_asym(:,2))
axis([0 200 0 100])
subplot(4,1,3)
plot(dist_asym(:,3))
axis([0 200 0 100])
subplot(4,1,4)
plot(dist_asym(:,4))
axis([0 200 0 100])