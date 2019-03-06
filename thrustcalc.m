gamma = 1.4 ;
R = 287.058 ; %[J/kg*K]
T_0 = 200 ; %[K]
h_p = 10000 ; %[m]
bpr = 2.6 ;
M_0 = V_0*sqrt(gamma*R*T_0) ;
m_f1 = 0.177 ;  %[kg/s] - fuel flow
m_f2 = 0.177 ;  %[kg/s] - fuel flow
m = 35.245 ;    %[kg/s] - total mass flow
T_ISA = T_0 + (bpr*h_p) ; %[K] - ISA Temperature

% ---------------------------------------------------------------------

data = rand(10,5) ;     % some random data
data = round(data*1000)/1000 ;
dlmwrite('matlab.dat',data,'delimiter',' ') 
load matlab.dat   % load the file
