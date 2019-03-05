
Dc=c/V0;
Db=b/V0;


C1_sym= [-2*muc*Dc, 0,0,0 ; 0, Dc*(CZadot-2*muc),0,0;0,0,-Db,0;0,Cmadot*Dc, 0, -2*muc*KY2*Dc];
C2_sym=[CXu,CXa,CZ0,CXq; CZu, CZa, -CX0, CZq+2*muc; 0,0,0,1; Cmu, Cma, 0, Cmq];
C3_sym=[CXde;  CZde; 0; Cmde];

C1_asym=[(CYbdot-2*mub)*Db, 0,0,0;0,-0.5*Db,0,0;0,0,-4*mub*KX2*Db, 4*mub*KXZ*Db; Cnbdot*Db,0,4*mub*KXZ*Db, -4*mub*KZ2*Db];
C2_asym=[CYb, CL, CYp, CYr; 0,0,1,0; Clb,0,Clp,Clr;Cnb,0,Cnp,Cnr];
C3_asym=[CYda, CYdr; 0,0; Clda, Cldr; Cnda, Cndr];


%define state space A and B
A_sym= -inv(C1_sym)*C2_sym;
B_sym= -inv(C1_sym)*C3_sym;

A_asym= -inv(C1_asym)*C2_asym;
B_asym= -inv(C1_asym)*C3_asym;

eig(A_sym)
eig(A_asym)


