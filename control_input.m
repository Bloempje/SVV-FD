% step iput 


figure
subplot(4,1,1)
plot(step(system_phugoid)(:,1))
subplot(4,1,2)
plot(step(system_phugoid)(:,2))
subplot(4,1,3)
plot(step(system_phugoid)(:,3))
subplot(4,1,4)
plot(step(system_phugoid)(:,4))
%legend('ú' 'a'  'thet' 'q')