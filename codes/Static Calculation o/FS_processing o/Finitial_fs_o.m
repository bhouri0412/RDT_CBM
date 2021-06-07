function Fi=Finitial_fs_o(Ug,Er,Gr,Iyy,Izz,Jt,alphap,Rr,Telt,Nbe,Npe,Kapfs)
% the effect of large displacement of ring freeshape on strain energy
% Npe=1000;

%element length
Le=Rr*Telt;

%discretization inside the element
dx=1/Npe;
x=0:dx:1;



%Shape functions and their derivatives
%========= for y and z displacements ===========%
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
%=========================================%

%========= for twist displacements ===========%
Nt1=1-3*x.^2+2*x.^3;
Nt2=Telt*(x-2*x.^2+x.^3);
Nt3=3*x.^2-2*x.^3;
Nt4=Telt*(-x.^2+x.^3);

dNt1=-6*x+6*x.^2;
dNt2=Telt*(1-4*x+3*x.^2);
dNt3=6*x-6*x.^2;
dNt4=Telt*(-2*x+3*x.^2);

Nt=[Nt1;Nt2;Nt3;Nt4];
dNt=[dNt1;dNt2;dNt3;dNt4];
%=========================================%
Fi=zeros(8*(Nbe+1),1);
Fiy=zeros(Nbe*Npe+1,1);
Fiz=zeros(Nbe*Npe+1,1);
Fit=zeros(Nbe*Npe+1,1);
for j=1:Nbe
    Fie=zeros(16,1);
    ind=1+(j-1)*Npe:1+j*Npe;
    indg=(8*(j-1)+1):(8*(j+1));
    
    %nodal y and z displacements for element j
    Uey=Ug([8*(j-1)+1:8*(j-1)+3,8*j+1:8*j+3]);
    Uez=Ug([8*(j-1)+4:8*(j-1)+6,8*j+4:8*j+6]);
    %nodal twist for element j
    Uet=Ug([8*j-1:8*j,8*j+7:8*j+8]);
    
    %ring radial displacement for element j
    yr=(Uey')*N;
    d2yr=(Uey')*d2N/Telt^2;
    %ring axial displacement for element j
    zr=(Uez')*N;
    dzr=(Uez')*dN/Telt;
    d2zr=(Uez')*d2N/Telt^2;
    %ring twist for element j
    alpha=(Uet')*Nt;
    dalpha=(Uet')*dNt/Telt;
    
    Kapyye=(1/Rr-(yr+d2yr)/Rr^2).*cos(alpha+d2zr/Rr);
    Kapzze=(1/Rr-(yr+d2yr)/Rr^2).*sin(alpha+d2zr/Rr);
    
    Kapfsyye=Kapfs(ind)'*cos(alphap);
    Kapfszze=Kapfs(ind)'*sin(alphap);
    
    dFiye=Le*Er*Izz*repmat((Kapyye-Kapfsyye).*(-cos(alpha+d2zr/Rr)/Rr^2),6,1).*(N+d2N/Telt^2)*dx...
        +Le*Er*Iyy*repmat((Kapzze-Kapfszze).*(-sin(alpha+d2zr/Rr)/Rr^2),6,1).*(N+d2N/Telt^2)*dx;
    Fiye=sum((dFiye(:,1:end-1)+dFiye(:,2:end))/2,2);
%     Dgye=sum(dDgye,2);
    
    dFize=Le*Er*Izz*repmat((Kapyye-Kapfsyye).*(-(1/Rr-(yr+d2yr)/Rr^2).*sin(alpha+d2zr/Rr)),6,1).*d2N/Telt^2/Rr*dx...
        +Le*Er*Iyy*repmat((Kapzze-Kapfszze).*(1/Rr-(yr+d2yr)/Rr^2).*cos(alpha+d2zr/Rr),6,1).*d2N/Telt^2/Rr*dx...
        +Le*Gr*Jt*repmat(dzr/Rr^2-dalpha/Rr,6,1).*dN/Telt/Rr^2*dx;
    Fize=sum((dFize(:,1:end-1)+dFize(:,2:end))/2,2);
%     Dgze=sum(dDgze,2);
    
    dFite=Le*Er*Izz*repmat((Kapyye-Kapfsyye).*(-(1/Rr-(yr+d2yr)/Rr^2).*sin(alpha+d2zr/Rr)),4,1).*Nt*dx...
        +Le*Er*Iyy*repmat((Kapzze-Kapfszze).*(1/Rr-(yr+d2yr)/Rr^2).*cos(alpha+d2zr/Rr),4,1).*Nt*dx...
        +Le*Gr*Jt*repmat(dzr/Rr^2-dalpha/Rr,4,1).*(-dNt/Telt/Rr)*dx;
    Fite=sum((dFite(:,1:end-1)+dFite(:,2:end))/2,2);
%     Dgte=sum(dDgte,2);
    
    Fie([1:3,9:11])=Fiye;
    Fie([4:6,12:14])=Fize;
    Fie([7,8,15,16])=Fite;
    Fi(indg)=Fi(indg)+Fie; % in N
end

% the=linspace(0,360,Nbe+1);
% 
% figure
% subplot(311)
% plot(the,Fi(1:8:end)/Le,'LineWidth',1.3)
% ylabel('dU/du_y (N/m)')
% set(gca,'xTick',[0:90:360])
% xlim([0,360])
% grid
% subplot(312)
% plot(the,Fi(4:8:end)/Le,'LineWidth',1.3)
% ylabel('dU/du_z (N/m)')
% set(gca,'xTick',[0:90:360])
% xlim([0,360])
% grid
% subplot(313)
% plot(the,Fi(7:8:end)/Le,'LineWidth',1.3)
% ylabel('dU/du_a (N)')
% xlabel('Circumferential (degree)')
% set(gca,'xTick',[0:90:360])
% xlim([0,360])
% grid
end
