function [Kap,r]=Kfs_cal_vg(Trfs,rfsi,rnod,Nbe,Npe,Telt)
%Calculate the interpolation of a given radius distribution given the
%nodal values and discretization. The curvature is also calculated

%discretization inside the element
dx=1/Npe;
x=(0:dx:1)';
Nr=Nbe*Npe;

%Shape functions
N1=1-10*x.^3+15*x.^4-6*x.^5;
N2=Telt*(x-6*x.^3+8*x.^4-3*x.^5);
N3=Telt^2*((x.^2)/2-3*(x.^3)/2+3*(x.^4)/2-(x.^5)/2);
N4=10*x.^3-15*x.^4+6*x.^5;
N5=Telt*(-4*x.^3+7*x.^4-3*x.^5);
N6=Telt^2*((x.^3)/2-x.^4+x.^5/2);

%first derivative of shape functions
dN1=-30*x.^2+60*x.^3-30*x.^4;
dN2=Telt*(1-18*x.^2+32*x.^3-15*x.^4);
dN3=Telt^2*(x-9*(x.^2)/2+6*x.^3-5*(x.^4)/2);
dN4=30*x.^2-60*x.^3+30*x.^4;
dN5=Telt*(-12*x.^2+28*x.^3-15*x.^4);
dN6=Telt^2*(3*(x.^2)/2-4*x.^3+5*(x.^4)/2);

%second derivative of shape functions
d2N1=-60*x+180*x.^2-120*x.^3;
d2N2=Telt*(-36*x+96*x.^2-60*x.^3);
d2N3=Telt^2*(1-9*x+18*x.^2-10*x.^3);
d2N4=60*x-180*x.^2+120*x.^3;
d2N5=Telt*(-24*x+84*x.^2-60*x.^3);
d2N6=Telt^2*(3*x-12*x.^2+10*x.^3);

%third derivative of shape functions
d3N1=-60+360*x-360*x.^2;
d3N2=Telt*(-36+192*x-180*x.^2);
d3N3=Telt^2*(-9+36*x-30*x.^2);
d3N4=60-360*x+360*x.^2;
d3N5=Telt*(-24+168*x-180*x.^2);
d3N6=Telt^2*(3-24*x+30*x.^2);

%intialization
r=zeros(Nr+1,1);
dr=zeros(Nr+1,1);
d2r=zeros(Nr+1,1);

for j=1:Nbe
rnode=rnod(3*(j-1)+1:3*j+3);
int=1+(j-1)*Npe:1+j*Npe;

%spline fonction
re=[N1,N2,N3,N4,N5,N6]*rnode;
%first derivative
dre=[dN1,dN2,dN3,dN4,dN5,dN6]*rnode/Telt;
%second derivative
d2re=[d2N1,d2N2,d2N3,d2N4,d2N5,d2N6]*rnode/Telt^2;
%third derivative
d3re=[d3N1,d3N2,d3N3,d3N4,d3N5,d3N6]*rnode/Telt^3;

%filling global matrices
r(int)=re;
dr(int)=dre;
d2r(int)=d2re;

end

%Calculation of curvature

Kap=(r.^2+2*dr.^2-r.*d2r)./(r.^2+dr.^2).^(3/2);

end

