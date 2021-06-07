function rdt_from_pr_to_fs_test_g()
close all
format long

%this code computes the free shape in millimeters from radial linear force as input: 
%the results are given within the variable named res_fs
%1st column contains angles in degrees and second one contains the radial coordinate of
%free shape in mm

%the code also gives the freeshape cuvature in the variable named res_Kapfs
%1st column contains angles in degrees and second one contains the free shape cuvature
%in 1/m

%for the input, 1st column should contain angles in degrees and second one force per
%circumference length in N/m 

Db=95.25e-3; %Bore diameter, m
Er=100e9; %Ring Young's modulus, Pa
Izr=10.6667e-12; %principal moment of inertia Izp, m^4
gap=0.314e-3; %ring gap size when closed inside bore
arm=2e-3; %running face minimum point width
Rr=Db/2-arm; %Ring nominal radius

% ask for text file containing two columns
%first column containing angles in degrees
%second column containing pressure times ring height (axial width) i.e.
%2.Lt/Dn*NormalizedPr
PressInput=dlmread('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/theo_pressure_TC_deg.txt');
Tpo=PressInput(:,1);
Pmo=PressInput(:,2);

Nbe=16; %number of elements
Npe=1000; %number of points within one element

Nr=Nbe*Npe;
Le=(2*pi*Rr-gap)/Nbe;
ds=Le/Npe;
Tgg=gap/2/Rr;
Tp=linspace(Tgg,pi,Nr/2+1)';
Tpend=linspace(pi,2*pi-Tgg,Nr/2+1)';
Pm=(interp1(Tpo,Pmo,Tp/pi*180,'linear','extrap')+flipud(interp1(Tpo,Pmo,Tpend/pi*180,'linear','extrap')))/2;

M=Ring_moment(Tp,Pm,Tgg,Rr);
kapfs=1/Rr-M/Er/Izr;
Kapfs1=[kapfs(1:end-1);flipud(kapfs)];
Tr=linspace(Tgg,2*pi-Tgg,Nbe*Npe+1)';

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
        M_ode=ringM(s,Pm,Tp,Tgg,Rr);
        
        feq(1)=yp(1)-sqrt(1-yp(2)^2)/y(2);
        feq(2)=yp(2)-y(3)*yp(1);
        feq(3)=1/Rr-M_ode/Er/Izr-(y(2)^2+2*y(3)^2-y(2)*yp(3)*sqrt(y(2)^2+y(3)^2))/(y(2)^2+y(3)^2)^(3/2);
        
    end

res_Kapfs=[Thefs*180/pi Kapfs1]; %Result giving the freeshape curvature
res_fs=[Thefs*180/pi Rfs]; %Result giving the free shape radial coordinates
dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/curvature_from_pressure_U_gap.txt',res_Kapfs,'delimiter','\t','precision','%12.8f','newline','pc');
dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/fs_from_pressure_U_gap.txt',res_fs,'delimiter','\t','precision','%12.8f','newline','pc');

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

function M_ode=ringM(s,Pm,Tp,Tg,Rr)

The=s/Rr;
alp=linspace(Tg,The,100)';
Pme=interp1(Tp,Pm,alp);
da=alp(2)-alp(1);

dM_ode=Pme*Rr^2.*sin(The-alp)*da;
M_ode=sum((dM_ode(1:end-1)+dM_ode(2:end))/2);

end

function M=Ring_moment(Tp,Pm,Tg,Rr)
The=linspace(Tp(1),Tp(end),180)';

Me=zeros(length(The),1);
for i=2:length(The)
    n=1000;
    Tha=linspace(Tg,The(i),n)';
    da=(The(i)-Tg)/n;
    Pmi=interp1(Tp,Pm,Tha);
    dMe=Pmi*Rr^2.*sin(The(i)-Tha)*da;
    Me(i)=sum(dMe(1:end-1)+dMe(2:end))/2;
end

M=interp1(The,Me,Tp);

end