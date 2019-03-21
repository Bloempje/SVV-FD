
mass=Mass_t*0.453592

%get all the initial value data at the stationary flight right before
%starting eigenmotion
[a0_ph, th0_ph, V_ph,q0_ph hp0_ph, m_ph, Db_ph, Dc_ph, muc_ph, mub_ph,CL_ph, CZ0_ph, CX0_ph]= flight_variables_sym(t1_phugoid,time, AoA, pitch, true_V, pressure_alt,mass, rho0,S,b,c, g, R, lambda, Temp0,q);
[a0_sp, th0_sp, V_sp,q0_sp, hp0_sp, m_sp, Db_sp, Dc_sp, muc_sp, mub_sp,CL_sp, CZ0_sp, CX0_sp]= flight_variables_sym(t1_shortperiod,time, AoA, pitch, true_V, pressure_alt,mass, rho0,S,b,c, g, R, lambda, Temp0,q);
[B0_dr,roll_angle0_dr, th0_dr, V_dr,p0_dr,r0_dr, hp0_dr, m_dr, Db_dr, Dc_dr, muc_dr, mub_dr,CL_dr, CZ0_dr, CX0_dr,index_dr]= flight_variables_asym(t1_dutchroll,time,sideslip,roll_angle, pitch, true_V, pressure_alt,mass, rho0,S,b,c, g, R, lambda, Temp0,p,r)
[B0_spiral,roll_angle0_spiral, th0_spiral, V_spiral,p0_spiral,r0_spiral, hp0_spiral, m_spiral,Db_spiral, Dc_spiral, muc_spiral, mub_spiral,CL_spiral, CZ0_spiral, CX0_spiral,index_sp]= flight_variables_asym(t1_spiral,time,sideslip,roll_angle, pitch, true_V, pressure_alt,mass, rho0,S,b,c, g, R, lambda, Temp0,p,r)
[B0_roll,roll_angle0_roll, th0_roll, V_roll,p0_roll,r0_roll, hp0_roll, m_roll, Db_roll, Dc_roll, muc_roll, mub_roll, CL_roll, CZ0_roll, CX0_roll,index_roll]= flight_variables_asym(t1_aperroll,time,sideslip,roll_angle, pitch, true_V, pressure_alt,mass,rho0,S,b,c, g, R, lambda, Temp0,p,r)

%plugging in all the variables for each eigenmotion and creating state
%space system
[system_phugoid,eigenvalues_phugoid, A_ph]=statespace_symmetrical(muc_ph,Dc_ph, V_ph, CZadot, Cmadot,Cma, KY2, CXu, CXa, CZ0_ph, CXq, CZu, CZa, CX0_ph, CZq, Cmu, Cmq, CXde,  CZde, Cmde)

[system_shortperiod,eigenvalues_shortperiod, A_sp]=statespace_symmetrical(muc_sp,Dc_sp, V_sp, CZadot, Cmadot,Cma, KY2, CXu, CXa, CZ0_sp, CXq, CZu, CZa, CX0_sp, CZq, Cmu, Cmq, CXde,  CZde, Cmde) 

[system_dutchroll, eigenvalues_dutchroll]= statespace_asymmetrical(mub_dr, V_dr , Db_dr, KX2, KXZ,KZ2, Cnbdot,CYbdot,CYb, CL_dr, CYp, CYr, Clb, Clp, Clr, Cnb, Cnp, Cnr, CYda, CYdr, Clda, Cldr, Cnda, Cndr )

[system_roll, eigenvalues_roll]= statespace_asymmetrical(mub_roll, V_roll , Db_roll, KX2, KXZ,KZ2, Cnbdot,CYbdot,CYb, CL_roll, CYp, CYr, Clb, Clp, Clr, Cnb, Cnp, Cnr, CYda, CYdr, Clda, Cldr, Cnda, Cndr )

[system_spiral, eigenvalues_spiral]= statespace_asymmetrical(mub_spiral, V_spiral , Db_spiral, KX2, KXZ,KZ2, Cnbdot,CYbdot,CYb, CL_spiral, CYp, CYr, Clb, Clp, Clr, Cnb, Cnp, Cnr, CYda, CYdr, Clda, Cldr, Cnda, Cndr )


%get the period and time to half amplitude based on eigenvalues of periodic
%dynamic motions
[P_phugoid, Thalf_phugoid]=periodic_characteristics(eigenvalues_phugoid(3))
[P_shortperiod, Thalf_shortperiod]=periodic_characteristics(eigenvalues_shortperiod(1))
[P_dutchroll, Thalf_dutchroll]=periodic_characteristics(eigenvalues_dutchroll(2))







%definitions needed to create state space system and eigenvalues

function [a0, th0, V,q0, hp0, m, Db, Dc, muc, mub, CL, CZ0, CX0]= flight_variables_sym(timemotion,time, alpha, pitch, V_true, pressure_alt, mass,  rho0,S,b,c, g, R, lambda, Temp0,q)
function  index=select_index(timemotion,time)
index=find(time==timemotion)
end
index=select_index(timemotion,time);
a0=alpha(index);
th0=pitch(index);
hp0=pressure_alt(index);
V=V_true(index);
m=mass(index)
q0=q(index)

rho = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1))
muc= m/(rho*S*c);
mub= m/(rho*S*b);
Dc=c/V;
Db=b/V;
W=m*g;
CL = 2*W/(rho*V^2*S); 
CZ0= -W*cos(deg2rad(th0))/(0.5*rho*V^2*S);
CX0= W*sin(deg2rad(th0))/(0.5*rho*V^2*S);
end 
function [B0,roll_angle0, th0, V,p0,r0 hp0, m, Db, Dc ,muc, mub, CL, CZ0, CX0,index]= flight_variables_asym(timemotion,time,sideslip,roll_angle, pitch, V_true, pressure_alt,mass,  rho0,S,b,c, g, R, lambda, Temp0,p,r)
function  index=select_index(timemotion,time)
index=find(time==timemotion)
end
index=select_index(timemotion,time);
B0=sideslip(index);
roll_angle0=roll_angle(index);
th0=pitch(index);
hp0=pressure_alt(index);
V=V_true(index);
m=mass(index);
p0=p(index)
r0=r(index)
rho = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1))
muc= m/(rho*S*c);
mub= m/(rho*S*b);
Dc=c/V;
Db=b/V;
W=m*g; 
CL = 2*W/(rho*V^2*S); 
CZ0=- W*cos(deg2rad(th0))/(0.5*rho*V^2*S);
CX0= W*sin(deg2rad(th0))/(0.5*rho*V^2*S);
end 



function [system_sym, sym_eigenvalues, A_sym]=statespace_symmetrical(muc,Dc, V, CZadot, Cmadot, Cma,  KY2, CXu, CXa, CZ0, CXq, CZu, CZa, CX0, CZq, Cmu, Cmq, CXde,  CZde, Cmde)
C1= [-2*muc*Dc/V, 0,0,0 ; 0, Dc*(CZadot-2*muc),0,0;0,0,-Dc,0;0,Cmadot*Dc, 0, -2*muc*KY2*Dc^2];
C2=[CXu/V,CXa,CZ0,CXq*Dc; CZu/V, CZa, -CX0, (CZq+2*muc)*Dc; 0,0,0,1*Dc; Cmu/V, Cma, 0, Cmq*Dc];
C3=[CXde;  CZde; 0; Cmde];

A_sym= -inv(C1)*C2;
B_sym= -inv(C1)*C3;
C_sym=eye(4);
D_sym=[0;0;0;0];

system_sym=ss(A_sym,B_sym,C_sym,D_sym);
sym_eigenvalues= eig(A_sym);
end 

function [system_asym, asym_eigenvalues]= statespace_asymmetrical(mub, V , Db, KX2, KXZ,KZ2, Cnbdot,CYbdot,CYb, CL, CYp, CYr, Clb, Clp, Clr, Cnb, Cnp, Cnr, CYda, CYdr, Clda, Cldr, Cnda, Cndr )
C1=[(CYbdot-2*mub)*Db, 0,0,0;0,-0.5*Db,0,0;0,0,-4*mub*KX2*Db^2/2, 4*mub*KXZ*Db^2/2; Cnbdot*Db,0,4*mub*KXZ*Db^2/2, -4*mub*KZ2*Db^2/2];
C2=[CYb, CL, CYp*Db/2, (CYr-4*mub)*Db/2; 0,0,1*Db/2,0; Clb,0,Clp*Db/2,Clr*Db/2;Cnb,0,Cnp*Db/2,Cnr*Db/2];
C3=[CYda, CYdr; 0,0; Clda, Cldr; Cnda, Cndr];
A_asym= -inv(C1)*C2;
B_asym= -inv(C1)*C3;
C_asym=eye(4);
D_asym=[0,0;0,0;0,0;0,0];
system_asym=ss(A_asym,B_asym,C_asym,D_asym) ;
asym_eigenvalues= eig(A_asym);
end 

function [P,Thalf]= periodic_characteristics(eigenvalue)
P=2*pi/(imag(eigenvalue))
Thalf=log(0.5)/(real(eigenvalue))
end











