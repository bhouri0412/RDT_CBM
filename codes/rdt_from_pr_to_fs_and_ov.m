function [fs,rov]=rdt_from_pr_to_fs_and_ov(inputs)

%this code computes the ovality in usual representation ring fixed 
%reference: 1st column contains angles in degrees and second one the 
%coordintes in microns in usual representation. it also computes free-shape
%in usual representation: 1st column contains angles in degrees in ring 
%fixed reference and second one the coordintes in millimeters  
%as inputs it needs a text file: 1st column contains angles in degrees in
%ring fixed reference and second one force par circumference length in N/m.

addpath(genpath(pwd))

format long

prompt='Bore diameter (mm): ';
Db=input(prompt)*1e-3;
prompt='Young modulus (GPa): ';
Er=input(prompt)*1e9; 
Izr=inputs.ring.Izr; 
prompt='Gap size when ring is closed to circular shape (mm): '; 
gap=input(prompt)*1e-3;
prompt='Gap size when ring is closed to ovality (under constant pressure) (mm): ';
gap2=input(prompt)*1e-3;
arm=inputs.ring.arm;
Rr=Db/2-arm;

%ask for file containing force distribution giving circular shape in ring
%fixed reference: first column containins angles in degrees and second one
%contains pressure times ring height (axial width) i.e. force per unit
%circumferential length in N/m
PressInput=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\pressure_TC.txt');  %user has to provide this
Tpo=PressInput(:,1);
Pmo=PressInput(:,2);
prompt='Type 1 to consider the constant linear force in computing the ovality equal to the average of the force distribution in radial bore provided, 0 to provide its value: ';
IsPovAvg=input(prompt);
if IsPovAvg
    Pov=mean(Pmo);
else
    prompt='Constant linear force (N/m) to consider in computing the ovality (Pressure times the ring height (axial width))';
    Pov=input(prompt);
end

prompt='Number of elements for the FEM: ';
Nbe=input(prompt); 
prompt='Number of points within one element: ';
Npe=input(prompt); 

prompt='Tolerance for Newton-Raphson algorithm convergence (suggested value: 1e-6): ';
Newton_tol=input(prompt);
prompt='Maximum number of Newton-Raphson algorithm iterations (suggested value: 100): ';
ItMax=input(prompt);

Nr=Nbe*Npe; 
Le=(2*pi*Rr-gap)/Nbe; 
ds=Le/Npe;
Tg=gap/2/Rr;
Tp=linspace(Tg,pi,Nr/2+1)';
Pm=interp1(Tpo,Pmo,Tp/pi*180,'linear','extrap');

M=Ring_moment(Tp,Pm,Tg,Rr);
kapfs=1/Rr-M/Er/Izr;
Kapfs=[kapfs(1:end-1);flipud(kapfs)];

s_span=[pi*Rr:-ds:Tg*Rr]';  

y0(1)=pi;     
y0(2)=Rr;     
y0(3)=0;      
yp0(1)=1/Rr; 
yp0(2)=0;     
yp0(3)=M(end)/Er/Izr;      

options=odeset('RelTol',1e-9);
[~,y,~]=ode15i(@fseqn,s_span',y0',yp0',options);

thefs=y(:,1);
Thefs=[flipud(thefs);2*pi-thefs(2:end)];
rfs=y(:,2);
Rfs=[flipud(rfs);rfs(2:end)];
fs=[Thefs*180/pi Rfs*1e3];

    function [feq]=fseqn(s,y,yp)
        feq=zeros(3,1);
        M_ode=ringM(s,Pm,Tp,Tg,Rr);
        
        feq(1)=yp(1)-sqrt(1-yp(2)^2)/y(2);
        feq(2)=yp(2)-y(3)*yp(1);
        feq(3)=1/Rr-M_ode/Er/Izr-(y(2)^2+2*y(3)^2-y(2)*yp(3)*sqrt(y(2)^2+y(3)^2))/(y(2)^2+y(3)^2)^(3/2);
        
    end

figure(1)
clf;
hold on
plot(fs(:,1),fs(:,2),'Linewidth',1.5)
xlabel('Circumferential Direction (degree)')
ylabel('Free shape in usual representation (mm)','FontSize',12)
set(gca,'XTick',[0:90:360])
xlim([0 360])
hold off
grid on
set(gca,'FontSize',16);

Tgg2=gap2/2/Rr;
Telt2=(2*pi-2*Tgg2)/(Nbe);
Kapfs=interp1(linspace(Tg,2*pi-Tg,Nr+1)',Kapfs,linspace(Tgg2,2*pi-Tgg2,Nr+1)','linear','extrap');
Tgg=Tgg2;
Telt=Telt2;
Le=Rr*Telt;

Ug=zeros(3*(Nbe+1),1);

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
    
    figure(10)
    cla
    semilogy(norm_NLS_k)
    grid
    Err_NLS_k=norm(NLS_k)
    
end

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

rov=[Tr' ov'*1e6];

xrnom=Rr*cosd(Tr);
yrnom=Rr*sind(Tr);
isnotSatisf=1;
while isnotSatisf
    prompt='Provide the magnifying coefficient for the radial plot of the ovality (recommanded values in order of 10): ';
    k=input(prompt);
    xov=(Rr+k*ov).*cosd(Tr);
    yov=(Rr+k*ov).*sind(Tr);
    kstr=num2str(k);
    s1='Plot of ovality in polar coordinates in meters with a magnifying coefficient equal to';
    s = [s1,' ',kstr];
    figure(2)
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
    grid
    axis equal
    set(gca,'FontSize',16);
    prompt='Type 1 if you want to choose again the magnifying coefficient, 0 otherwise: ';
    isnotSatisf=input(prompt);
end

figure(3)
clf;
hold on
plot(Tr',ov'*1e6,'LineWidth',1.5)
xlabel('Circumferential Direction (degree)')
ylabel('Ovality (\mu m)')
title('Ovality in usual representation');
set(gca,'XTick',[0:90:360])
set(gca,'FontSize',16);
grid on

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