function [E_stac,yr,zr,alr,hmin,fl,z0,el1,el2,el3]=Static_Error_vg(Ug_nd,F1,inputs)
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
[Fext_nd,yr,zr,alr,hmin,fl,z0]=Ext_Forces_static_vg(Ug_nd,inputs);

F1_nd=zeros(8*Nbnod,1);
F1_nd(1:8:end-7)=F1(1:8:end-7)*Rr^4*Telt^4/Le/Er/Izz/uref;
F1_nd(2:8:end-6)=F1(2:8:end-6)*Rr^4*Telt^3/Le/Er/Izz/uref;
F1_nd(3:8:end-5)=F1(3:8:end-5)*Rr^4*Telt^2/Le/Er/Izz/uref;

F1_nd(4:8:end-4)=F1(4:8:end-4)*Rr^4*Telt^4/Le/Er/Iyy/href;
F1_nd(5:8:end-3)=F1(5:8:end-3)*Rr^4*Telt^3/Le/Er/Iyy/href;
F1_nd(6:8:end-2)=F1(6:8:end-2)*Rr^4*Telt^2/Le/Er/Iyy/href;

F1_nd(7:8:end-1)=F1(7:8:end-1)*Rr^3*Telt^2/Le/Gr/Jt/href;
F1_nd(8:8:end)=F1(8:8:end)*Rr^3*Telt/Le/Gr/Jt/href;

%Calculating the system residual
E_stac=Kg_nd*Ug_nd-(Fext_nd-F1_nd);
el1=Kg_nd*Ug_nd;
el2=Fext_nd;
el3=F1_nd;
end