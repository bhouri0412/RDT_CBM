function [Thefs,Rfs,Kapfs]=fscal_n4g(Tpo,Pmo,Rr,Er,Izr,gap,Nbe,Npe)

% FEM input
Nr=Nbe*Npe; %compatibility
Le=(2*pi*Rr-gap)/Nbe; %Element length (mm)
ds=Le/Npe;
Tgg=gap/2/Rr;
% Tg=0;
% save_step=round(Nbe*Npe/360);

Tp=linspace(Tgg,pi,Nr/2+1)';
Pm=interp1(Tpo,Pmo,Tp/pi*180);

M=Ring_moment(Tp,Pm,Tgg,Rr);
kapfs=1/Rr-M/Er/Izr;
Kapfs=[kapfs(1:end-1);flipud(kapfs)];
% Kapov=Kapfs+Mov/Er/Izr;

s_span=[pi*Rr:-ds:Tgg*Rr]';     %from 180 to Tg

y0(1)=pi;     %theta
y0(2)=Rr;     %r
y0(3)=0;      %dr/dth
yp0(1)=1/Rr;  %dth/ds
yp0(2)=0;     %dr/ds
yp0(3)=M(end)/Er/Izr/Rr;       %d(dr/dth)/ds

options=odeset('RelTol',1e-9);
[s,y,yp]=ode15i(@fseqn,s_span',y0',yp0',options);

thefs=y(:,1);
Thefs=[flipud(thefs);2*pi-thefs(2:end)];
rfs=y(:,2);
Rfs=[flipud(rfs);rfs(2:end)];
drfs=y(:,3);
dRfs=[flipud(drfs);-drfs(2:end)];
fs=[Thefs/pi*180 Rfs*1e6];

function [feq]=fseqn(s,y,yp)
    feq=zeros(3,1);
    M_ode=ringM(s,Pm,Tp,Tgg,Rr);
    
    feq(1)=yp(1)-sqrt(1-yp(2)^2)/y(2);
    feq(2)=yp(2)-y(3)*yp(1);
    feq(3)=1/Rr-M_ode/Er/Izr-(y(2)^2+2*y(3)^2-y(2)*yp(3)*sqrt(y(2)^2+y(3)^2))/(y(2)^2+y(3)^2)^(3/2);

end

end
function M_ode=ringM(s,Pm,Tp,Tg,Rr)

The=s/Rr;
alp=linspace(Tg,The,1000)';
Pme=interp1(Tp,Pm,alp);
da=alp(2)-alp(1);

dM_ode=Pme*Rr^2.*sin(The-alp)*da;
M_ode=sum((dM_ode(1:end-1)+dM_ode(2:end))/2);

end

function M=Ring_moment(Tp,Pm,Tg,Rr)
The=linspace(Tp(1),Tp(end),1800)';

Me=zeros(length(The),1);
for i=2:length(The)
    n=100000;
    Tha=linspace(Tg,The(i),n)';
    da=(The(i)-Tg)/n;
    Pmi=interp1(Tp,Pm,Tha);
    dMe=Pmi*Rr^2.*sin(The(i)-Tha)*da;
    Me(i)=sum(dMe(1:end-1)+dMe(2:end))/2;
end

M=interp1(The,Me,Tp);

end
function plot_freeshape(fs)
fs_rdt=dlmread('C:\Libraries\MIT Work\Models\Top2Ring oil transfer model\Static Calculation\test_results\rdt\fscal\uniform_32.txt');

figure,
hold on
plot(fs(:,1),fs(:,2)/1e3,'LineWidth',1.5)
% plot(fs(1:Npe:end,1)/pi*180,fs(1:Npe:end,2)*1e3,'ro')
plot(fs_rdt(:,1),fs_rdt(:,2)/1e3,'--','Color',[0 0.8 0],'LineWidth',1.5)
xlabel('Circumferential Direction (degree)')
ylabel('Radial Coordinate(mm)')
set(gca,'XTick',[0:90:360])
xlim([0,360])
grid
end