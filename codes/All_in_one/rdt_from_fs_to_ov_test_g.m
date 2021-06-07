function rdt_from_fs_to_ov_test_g()
close all
format long

%this code computes the ovality in usual representation: 
%the results are given within the variable named res_ov
%1st column contains angles in degrees and second one coordintes in microns

%the user either provides the free shape coordinates with first column
%containing angles in degrees and the second one containing the freeshape
%in millimeters
%or provides the free shape curvature with first column containing angles
%in degrees and the second one containing the freeshape curvature in 1/m

Db=95.25e-3;  %bore diameter, m
Er=100e9; %Young's modulus, Pa, ring
Izr=10.6667e-12; %principal moment of inertia Izp, m^4
arm=2e-3; %running face minimum point width
gap=0.314e-3; %ring gap size when closed to ovality under constant pressure
IsFreeshape=0; %1 if user provides free shape coordinates, 0 if provides free shape curvature

%the constant radial force needs to be provided in the variable Pov in N/m
%it is the pressure times ring height (axial width)
PressInput=dlmread('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/theo_pressure_U_deg.txt');
Pmo=PressInput(:,2);
Pov=mean(Pmo);

Rr=Db/2-arm;

Nbe=16; %number of elements
Npe=1000; %number of points within a element for contact calculation
Le=(2*pi*Rr-gap)/Nbe;
Tgg=gap/2/Rr;

if IsFreeshape %considering free shape coordinates

    %ask for free shape data with first column containing angles in degrees and
    %the second one containing the freeshape in millimeters

    fs=dlmread('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/fs_from_theo_pr_TC.txt');

    Tri=fs(:,1)*pi/180;
    rfsim=fs(:,2)*1e-3;
    Inod=1:Npe:Npe*Nbe+1;
    dTrfs=(Tri(end)-Tri(1))/(Npe*Nbe);
    Trfs=(Tri(1):dTrfs:Tri(length(Tri)))';
    Teltfs=(Trfs(end)-Trfs(1))/Nbe;
    Tnodfs=Trfs(Inod);
    rfsi=interp1(Tri,rfsim,Trfs);
    [Ccs]=Cs_3rd_cont_o(Tnodfs);
    [rnodfs]=ls_rfsi_o(rfsi,Nbnod,Teltfs,Ccs);
    [Kapfs,~]=Kfs_cal_o(Trfs,rfsi,rnodfs,Nbe,Npe,Teltfs);
else
    
    %ask for free shape curvature with first column containing angles in degrees and
    %the second one containing the freeshape curvature in 1/m
    Kapfs=dlmread('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/Kapfs_from_theo_pr_U_gap.txt');
    
    Tri=Kapfs(:,1)*pi/180;
    Kapfs=interp1(Tri,Kapfs(:,2),linspace(Tgg,2*pi-Tgg,Npe*Nbe+1)','linear','extrap');
end

Ug=zeros(3*(Nbe+1),1);

Newton_tol=1e-6; %Tolerance for Newton-Raphson algorithm convergence (suggested value: 1e-6)
ItMax=100; %Maximum number of Newton-Raphson algorithm iterations (suggested value: 100)

judge1=1;
NLS_0=NewtonOV_Error(Ug,Er,Izr,Rr,Nbe,Npe,Le,Kapfs',Pov);
NLS_k=NLS_0;
j=0;

while j<=ItMax&&(judge1)
    j=j+1;

    [Ja_ov]=Jacobian_ov(Nbe,Npe,Le,Er,Izr,Rr,Ug,Kapfs');
    Nstep=-Ja_ov\NLS_k;
    Ug=Ug+Nstep;

    NLS_k=NewtonOV_Error(Ug,Er,Izr,Rr,Nbe,Npe,Le,Kapfs',Pov);
    judge1=norm(NLS_k)>Newton_tol;
    norm_NLS_k(j)=norm(NLS_k);
    
    figure(1)
    cla
    semilogy(norm_NLS_k)
    grid
    pause(0.0001)
    Err_NLS_k=norm(NLS_k)
end

res_ov=plot_ovality(Ug,Rr,Le,Nbe,Npe,Tgg);

end

function res_ov=plot_ovality(Ug,Rr,Le,Nbe,Npe,Tgg)
dx=1/Npe;
x=0:dx:1;
Tr=linspace(Tgg*180/pi,360-Tgg*180/pi,Nbe*Npe+1);

Telt=Le/Rr;

N1=1-10*x.^3+15*x.^4-6*x.^5;
N2=Telt*(x-6*x.^3+8*x.^4-3*x.^5);
N3=Telt^2*((x.^2)/2-3*(x.^3)/2+3*(x.^4)/2-(x.^5)/2);
N4=10*x.^3-15*x.^4+6*x.^5;
N5=Telt*(-4*x.^3+7*x.^4-3*x.^5);
N6=Telt^2*((x.^3)/2-x.^4+x.^5/2);
N=[N1;N2;N3;N4;N5;N6];

ov=zeros(1,Nbe*Npe+1);
for j=1:Nbe
    ind=1+(j-1)*Npe:1+j*Npe;
    indg=(3*(j-1)+1):(3*(j+1));
    
    Uey=Ug(indg);
    yr=(Uey')*N;
    ov(ind)=yr;
end

res_ov=[Tr' 1e6*ov'];
dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/Theo_OV_U_gap.txt',res_ov,'delimiter','\t','precision','%12.8f','newline','pc');

xrnom=Rr*cosd(Tr);
yrnom=Rr*sind(Tr);
k=30; %magnifying coefficient
xov=(Rr+k*ov).*cosd(Tr);
yov=(Rr+k*ov).*sind(Tr);

kstr=num2str(k);
s1='Ovality plot with a magnifying coefficient equal to';
s = [s1,' ',kstr];
figure(1)
clf;
hold on
plot((xov-(max(xov)+min(xov))/2)*1e3,(yov-(max(yov)+min(yov))/2)*1e3,'-r','Linewidth',1.5)
plot(xrnom*1e3,yrnom*1e3,'--k')
plot((xov(1:Npe:end)-(max(xov)+min(xov))/2)*1e3,(yov(1:Npe:end)-(max(yov)+min(yov))/2)*1e3,'or','Linewidth',1.5)
legend('Ovality','Ring nominal radius','Location','northoutside')
title(s);
xlabel('mm');
ylabel('mm');
hold off
grid on
axis equal
set(gca,'FontSize',16);

figure(2)
clf;
hold on
plot(Tr,ov*1e6,'LineWidth',1.5)
ylabel('Ovality (\mu m)')
set(gca,'XTick',[0:90:360])
set(gca,'FontSize',16);
grid on

end

function NLS_k=NewtonOV_Error(Ug,Er,Izr,Rr,Nbe,Npe,Le,Kapfs,Pov)
dx=1/Npe;
x=0:dx:1;

Telt=Le/Rr;

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
Fm=zeros(3*(Nbe+1),1);
dU=zeros(3*(Nbe+1),1);

for j=1:Nbe
    ind=1+(j-1)*Npe:1+j*Npe;
    indg=(3*(j-1)+1):(3*(j+1));
    
    Uey=Ug(indg);
    
    yr=(Uey')*N;
    dyr=(Uey')*dN/Telt;
    d2yr=(Uey')*d2N/Telt^2;
   
    r=Rr+yr;

    kap=(r.^2-r.*d2yr+2*dyr.^2)./(r.^2+dyr.^2).^(3/2);
    dkapdy=(2*r-d2yr)./(r.^2+dyr.^2).^(3/2)-3*(r.^3-r.^2.*d2yr+2*r.*dyr.^2)./(r.^2+dyr.^2).^(5/2);
    dkapddy=4*dyr./(r.^2+dyr.^2).^(3/2)-3*(r.^2-r.*d2yr+2*dyr.^2).*dyr./(r.^2+dyr.^2).^(5/2);
    dkapdd2y=-r./(r.^2+dyr.^2).^(3/2);

    dkap=kap-Kapfs(ind);
    dkapn=repmat(dkap,6,1);
    
    dkapdyn=repmat(dkapdy,6,1);
    dkapddyn=repmat(dkapddy,6,1);
    dkapdd2yn=repmat(dkapdd2y,6,1);
    
    ddUe=Le*Er*Izr*dkapn.*(dkapdyn.*N+dkapddyn.*dN/Telt+dkapdd2yn.*d2N/Telt^2)*dx;
    dUe=sum((ddUe(:,1:end-1)+ddUe(:,2:end))/2,2);
    
    dFme=-Le*Pov*N*dx;
    Fme=sum((dFme(:,1:end-1)+dFme(:,2:end))/2,2);
    
 
    dU(indg)=dU(indg)+dUe;
    Fm(indg)=Fm(indg)+Fme;
end

NLS_k=dU-Fm;
end

function [Ja_ov]=Jacobian_ov(Nbe,Npe,Le,Er,Izr,Rr,Ug,Kapfs)

Telt=Le/Rr;

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

Ja_ov=zeros(3*(Nbe+1),3*(Nbe+1));

for n=1:Nbe
    Ja_ove=zeros(6,6);
    ind=1+(n-1)*Npe:1+n*Npe;
    indg=(3*(n-1)+1):(3*(n+1));
    
    Uey=Ug(indg);
    
    yr=(Uey')*N;
    dyr=(Uey')*dN/Telt;
    d2yr=(Uey')*d2N/Telt^2;

    r=Rr+yr;

    kap=(r.^2-r.*d2yr+2*dyr.^2)./(r.^2+dyr.^2).^(3/2);
    
    dkapdy=(2*r-d2yr)./(r.^2+dyr.^2).^(3/2)-3*(r.^3-r.^2.*d2yr+2*r.*dyr.^2)./(r.^2+dyr.^2).^(5/2);
    dkapddy=4*dyr./(r.^2+dyr.^2).^(3/2)-3*(r.^2-r.*d2yr+2*dyr.^2).*dyr./(r.^2+dyr.^2).^(5/2);
    dkapdd2y=-r./(r.^2+dyr.^2).^(3/2);
    
    d2kapdydy=2./(r.^2+dyr.^2).^(3/2)-(15*r.^2-9*d2yr.*r+6*dyr.^2)./(r.^2+dyr.^2).^(5/2)...
        +15*(r.^4-r.^3.*d2yr+2*r.^2.*dyr.^2)./(r.^2+dyr.^2).^(7/2);
    d2kapdyddy=(-18*r.*dyr+3*dyr.*d2yr)./(r.^2+dyr.^2).^(5/2)...
        +15*(r.^3.*dyr-r.^2.*dyr.*d2yr+2*r.*dyr.^3)./(r.^2+dyr.^2).^(7/2);
    d2kapdydd2y=-(r.^2+dyr.^2).^(-3/2)+3*r.^2./(r.^2+dyr.^2).^(5/2);
    
    d2kapddyddy=4*(r.^2+dyr.^2).^(-3/2)-(30*dyr.^2+3*r.^2-3*r.*d2yr)./(r.^2+dyr.^2).^(5/2)...
        +15*(r.^2.*dyr.^2-r.*dyr.^2.*d2yr+2*dyr.^4)./(r.^2+dyr.^2).^(7/2);
    d2kapddydd2y=3*r.*dyr./(r.^2+dyr.^2).^(5/2);

    dkap=kap-Kapfs(ind);
    
    for i=1:6
        for j=1:6
            
            dJa_fsije=Le*Er*Izr*((dkapdy.*N(j,:)+dkapddy.*dN(j,:)/Telt+dkapdd2y.*d2N(j,:)./Telt^2).*(dkapdy.*N(i,:)+dkapddy.*dN(i,:)/Telt+dkapdd2y.*d2N(i,:)./Telt^2)...
                +dkap.*(d2kapdydy.*N(j,:).*N(i,:)+d2kapdyddy.*dN(j,:)/Telt.*N(i,:)+d2kapdydd2y.*d2N(j,:).*N(i,:)/Telt^2 ...
                +d2kapdyddy.*N(j,:).*dN(i,:)/Telt+d2kapddyddy.*dN(j,:).*dN(i,:)/Telt^2+d2kapddydd2y.*d2N(j,:).*dN(i,:)/Telt^3 ...
                +d2kapdydd2y.*N(j,:).*d2N(i,:)/Telt^2+d2kapddydd2y.*dN(j,:).*d2N(i,:)/Telt^3))*dx;
            
            Ja_fsije1=sum((dJa_fsije(:,1:end-1)+dJa_fsije(:,2:end))/2,2);
            
           
            
            Ja_ove(i,j)=Ja_fsije1;

            
        end
    end
    Ja_ov(indg,indg)=Ja_ov(indg,indg)+Ja_ove;

end
end
