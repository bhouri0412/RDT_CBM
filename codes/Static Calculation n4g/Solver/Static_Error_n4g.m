function [E_stac,yr,zr,alr,hmin,fl,z0,el1,el2,el3]=Static_Error_n4g(Ug_nd,F1,inputs)
%Error function of the dynamics equation of the ring FEM
Rr=inputs.ring.Rr; %ring neutral fiber radius
Le=inputs.ring.Le;
Nbnod=inputs.FEM.Nbnod;
Telt=inputs.FEM.Telt; %element size (rad)
uref=inputs.FEM.uref;
href=inputs.FEM.href;
Er=inputs.ring.Er;
Gr=inputs.ring.Gr;
Izz=inputs.ring.Izz;
Iyy=inputs.ring.Iyy;
Jt=inputs.ring.Jt;

%Retrieving FEM matrices from input structure
Kg_nd=inputs.FEM.Kg_nd;

%calculation of external nodal forces
[Fext_nd,yr,zr,alr,hmin,fl,z0]=Ext_Forces_static_n4g(Ug_nd,inputs);

F1_nd=zeros(12*Nbnod,1);
F1_nd(1:12:end-11)=F1(1:12:end-11)*Rr^4*Telt^4/Le/Er/Izz/uref;
F1_nd(2:12:end-10)=F1(2:12:end-10)*Rr^4*Telt^3/Le/Er/Izz/uref;
F1_nd(3:12:end-9)=F1(3:12:end-9)*Rr^4*Telt^2/Le/Er/Izz/uref;

F1_nd(4:12:end-8)=F1(4:12:end-8)*Rr^4*Telt^4/Le/Er/Iyy/href;
F1_nd(5:12:end-7)=F1(5:12:end-7)*Rr^4*Telt^3/Le/Er/Iyy/href;
F1_nd(6:12:end-6)=F1(6:12:end-6)*Rr^4*Telt^2/Le/Er/Iyy/href;

F1_nd(7:12:end-5)=F1(7:12:end-5)*Rr^3*Telt^2/Le/Gr/Jt/href;
F1_nd(8:12:end-4)=F1(8:12:end-4)*Rr^3*Telt/Le/Gr/Jt/href;
F1_nd(9:12:end-3)=F1(9:12:end-3)*Rr^3/Le/Gr/Jt/href;

%Calculating the system residual
E_stac=Kg_nd*Ug_nd-(Fext_nd-F1_nd);
el1=Kg_nd*Ug_nd;
el2=Fext_nd;
el3=F1_nd;
end