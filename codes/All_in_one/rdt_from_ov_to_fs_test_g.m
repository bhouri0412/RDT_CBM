function rdt_from_ov_to_fs_test_g()

%this code computes the free shape from ovality in usual representation as input: 
%the results are given within the variable named res_fs
%1st column contains angles in degrees and second one contains the radial coordinate of
%free shape in mm

%the code also gives the freeshape cuvature in the variable named res_Kapfs
%1st column contains angles in degrees and second one contains the free shape cuvature
%in 1/m

%for the input, 1st column should contain angles in degrees and second one ovality radial
%coordinates in microns

close all

Db=95.25e-3;  %bore diameter, m
Er=100e9; %Young's modulus, Pa, ring
Izr=10.6667e-12; %principal moment of inertia Izp, m^4
arm=2e-3; %running face minimum point width
gap=0.314e-3; %ring gap size when closed to ovality under constant pressure

Nbe=16; %number of elements
Npe=1000; %number of points within a element for contact calculation

%ask for ovality input with angles in degrees in first column and ovality
%in usual representation in microns
ovdata=dlmread('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/Theo_OV_U_gap.txt');

%the constant radial force needs to be provided in the variable Pov in N/m
%it is the pressure times ring height (axial width)
PressInput=dlmread('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/theo_pressure_U_deg.txt');
Pmo=PressInput(:,2);
Pov=mean(Pmo);

Rr=Db/2-arm;
Tgg=gap/2/Rr;

Telt=(2*pi-2*Tgg)/(Nbe);
Tr=linspace(Tgg,2*pi-Tgg,Nbe*Npe+1)';

ovdata(:,1)=ovdata(:,1)*pi/180;
ovdata(:,2)=ovdata(:,2)*1e-6;
Uov=ov_lsm_check_o(ovdata,Nbe,Npe);

ov=zeros(1,Npe*Nbe+1);
dov=zeros(1,Npe*Nbe+1);
d2ov=zeros(1,Npe*Nbe+1);

dx=1/Npe;
x=0:dx:1;

N1=1-10*x.^3+15*x.^4-6*x.^5;
N2=Telt*(x-6*x.^3+8*x.^4-3*x.^5);
N3=Telt^2*((x.^2)/2-3*(x.^3)/2+3*(x.^4)/2-(x.^5)/2);
N4=10*x.^3-15*x.^4+6*x.^5;
N5=Telt*(-4*x.^3+7*x.^4-3*x.^5);
N6=Telt^2*((x.^3)/2-x.^4+x.^5/2);

dN1=-30*x.^2+60*x.^3-30*x.^4;
dN2=Telt*(1-18*x.^2+32*x.^3-15*x.^4);
dN3=Telt^2*(x-9*(x.^2)/2+6*x.^3-5*(x.^4)/2);
dN4=30*x.^2-60*x.^3+30*x.^4;
dN5=Telt*(-12*x.^2+28*x.^3-15*x.^4);
dN6=Telt^2*(3*(x.^2)/2-4*x.^3+5*x.^4/2);

d2N1=-60*x+180*x.^2-120*x.^3;
d2N2=Telt*(-36*x+96*x.^2-60*x.^3);
d2N3=Telt^2*(1-9*x+18*x.^2-10*x.^3);
d2N4=60*x-180*x.^2+120*x.^3;
d2N5=Telt*(-24*x+84*x.^2-60*x.^3);
d2N6=Telt^2*(3*x-12*x.^2+10*x.^3);

N=[N1;N2;N3;N4;N5;N6];
dN=[dN1;dN2;dN3;dN4;dN5;dN6];
d2N=[d2N1;d2N2;d2N3;d2N4;d2N5;d2N6];

for i=1:Nbe-1
    
    Uey=Uov(3*i-2:3*i+3,1);
    
    ov(1+(i-1)*Npe:i*Npe)=(Uey')*N(:,1:end-1);
    dov(1+(i-1)*Npe:i*Npe)=(Uey')*dN(:,1:end-1)/Telt;
    d2ov(1+(i-1)*Npe:i*Npe)=(Uey')*d2N(:,1:end-1)/Telt^2;
end

Uey=Uov(end-5:end);
ov(1+(Nbe-1)*Npe:1+Nbe*Npe)=(Uey')*N;
dov(1+(Nbe-1)*Npe:1+Nbe*Npe)=(Uey')*dN/Telt;
d2ov(1+(Nbe-1)*Npe:1+Nbe*Npe)=(Uey')*d2N/Telt^2;

ov=ov';
dov=dov';
d2ov=d2ov';

Nr=Nbe*Npe; 
Le=(2*pi*Rr-gap)/Nbe;
ds=Le/Npe;
kov=((Rr+ov(1:Nr/2+1)).^2+2*dov(1:Nr/2+1).^2-(Rr+ov(1:Nr/2+1)).*d2ov(1:Nr/2+1))./((Rr+ov(1:Nr/2+1)).^2+dov(1:Nr/2+1).^2).^(3/2);
Tp=linspace(Tgg,pi,Nr/2+1)'; 
Pm=Pov*ones(length(Tp),1);

M=Ring_moment(Tp,Pm,Tr,ov,dov,Tgg,Rr);
kapfs=kov-M/Er/Izr;
Kapfs1=[kapfs(1:end-1);flipud(kapfs)];

s_span=[pi*Rr:-ds:Tgg*Rr]';

y0(1)=pi; 
y0(2)=Rr;  
y0(3)=0;     
yp0(1)=1/Rr; 
yp0(2)=0;    
yp0(3)=M(end)*Rr/Er/Izr;     

options=odeset('RelTol',1e-9);
[~,y,~]=ode15i(@fseqn,s_span',y0',yp0',options);

thefs=y(:,1);
Thefs=[flipud(thefs);2*pi-thefs(2:end)];
rfs=y(:,2);
Rfs=1e3*[flipud(rfs);rfs(2:end)];

function [feq]=fseqn(s,y,yp)
    feq=zeros(3,1);
    M_ode=ringM(s,Pm,Tp,Tgg,ov,dov,Rr,Tr,Npe,Nbe);
    K_ode=ringK(s,ov,dov,d2ov,Rr,Npe,Nbe);
    feq(1)=yp(1)-sqrt(1-yp(2)^2)/y(2);
    feq(2)=yp(2)-y(3)*yp(1);
    feq(3)=K_ode-M_ode/Er/Izr-(y(2)^2+2*y(3)^2-y(2)*yp(3)*sqrt(y(2)^2+y(3)^2))/(y(2)^2+y(3)^2)^(3/2);

end

res_Kapfs=[Thefs*180/pi Kapfs1]; %Result giving the freeshape curvature
res_fs=[Thefs*180/pi Rfs]; %Result giving the free shape radial coordinates
dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/curvature_from_OV_U_gap.txt',res_Kapfs,'delimiter','\t','precision','%12.8f','newline','pc');
dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/fs_from_OV_U_gap.txt',res_fs,'delimiter','\t','precision','%12.8f','newline','pc');

figure(1)
clf;
hold on
plot(Thefs*180/pi,Kapfs1,'LineWidth',1.5)
xlabel('Circumferential Direction (degree)')
ylabel('Curvature of free shape (1/m)')
xlim([0 360])
set(gca,'XTick',[0:90:360])
grid on
set(gca,'FontSize',16);

figure(2)
clf;
hold on
plot(Thefs*180/pi,Rfs,'LineWidth',1.5)
xlabel('Circumferential Direction (degree)')
ylabel('Free shape (mm)')
xlim([0 360])
set(gca,'XTick',[0:90:360])
grid on
set(gca,'FontSize',16);

xrnom=1e3*Rr*cos(linspace(0,2*pi,Nbe*Npe+1)');
yrnom=1e3*Rr*sin(linspace(0,2*pi,Nbe*Npe+1)');

xfsc=Rfs.*cos(Thefs);
yfsc=Rfs.*sin(Thefs);

figure(3)
clf;
hold on
plot(xfsc,yfsc,'LineWidth',1.5)
plot(xrnom,yrnom,'--k','LineWidth',1.5)
legend('Free shape','Ring nominal radius','Location','northoutside')
title('Free shape (in mm)')
axis equal
set(gca,'xtick',[-50:10:50])
set(gca,'ytick',[-50:10:50])
xlim([-50,50])
ylim([-50,50])
grid on
set(gca,'FontSize',16);

end

function K_ode=ringK(s,ov,dov,d2ov,Rr,Npe,Nbe)
i=floor((pi*Rr-s)*Npe*Nbe/(2*pi*Rr))+1;
j=Npe*Nbe/2+2-i;
K_ode=((Rr+ov(j,1))^2+2*dov(j,1)^2-(Rr+ov(j,1))*d2ov(j,1))/((Rr+ov(j,1))^2+dov(j,1)^2)^(3/2);
end
function M_ode=ringM(s,Pm,Tp,Tg,ov,dov,Rr,Tr,Npe,Nbe)

i=floor((pi*Rr-s)*Npe*Nbe/(2*pi*Rr))+1;
j=Npe*Nbe/2+2-i;
if j>2
The=Tr(j,1);
alp=linspace(Tg,The,1000)';

uov=interp1(Tr,ov,alp);
udov=interp1(Tr,dov,alp);
upov=interp1(Tp,Pm,alp);

df=abs((udov.*cos(alp)-(Rr+uov).*sin(alp)).*(Rr+ov(j,1)).*cos(The)+(udov.*sin(alp)+(Rr+uov).*cos(alp)).*(Rr+ov(j,1)).*sin(The)-(Rr+uov).*udov)./sqrt((udov.*cos(alp)-(Rr+uov).*sin(alp)).^2+(udov.*sin(alp)+(Rr+uov).*cos(alp)).^2);
da=alp(2)-alp(1);

dM_ode=upov.*(Rr+uov).*df*da;
M_ode=sum((dM_ode(1:end-1)+dM_ode(2:end))/2);
else
    M_ode=0;
end
end

function M=Ring_moment(Tp,Pm,Tr,ov,dov,Tg,Rr)

M=zeros(length(Tp),1);

for i=2:length(Tp)
    n=10000;
    Tha=linspace(Tg,Tp(i),n)';
    uov=interp1(Tr,ov,Tha);
    udov=interp1(Tr,dov,Tha);
    da=(Tp(i)-Tg)/n;
    upov=interp1(Tp,Pm,Tha);
    df=abs((udov.*cos(Tha)-(Rr+uov).*sin(Tha)).*(Rr+ov(i,1)).*cos(Tp(i))+(udov.*sin(Tha)+(Rr+uov).*cos(Tha)).*(Rr+ov(i,1)).*sin(Tp(i))-(Rr+uov).*udov)./sqrt((udov.*cos(Tha)-(Rr+uov).*sin(Tha)).^2+(udov.*sin(Tha)+(Rr+uov).*cos(Tha)).^2);
    dMe=upov.*(Rr+uov).*df*da;
    M(i)=sum(dMe(1:end-1)+dMe(2:end))/2;
end


end
