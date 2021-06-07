function rdt_conformability_radial_plot(Nplot,Nr,ybr,Db,Tg,hminv,fyring,flring,IsBDist,gap)

%this code generates the radial plots related to the conformability for one
%specific ring gap position once the calculations are carried out via the
%function rdt_conformability

addpath(genpath(pwd))
dtg=180*gap/pi/2/Rr+0.2;
prompt='Magnification coefficient for radial plot of bore distortion and ring liner clearance (suggested value: 3000): ';
mag1=input(prompt);
mag2=mag1;
prompt='Number of radial forces to represent (suggested value: 200): ';
numforce=input(prompt);
prompt='Maximum acceptable value for forces to be represented (N/m) (suggested value: 20000): ';
fstop=input(prompt); 

if IsBDist
    
theeplotv=linspace(0,360,Nplot+1);
theev=linspace(0,360,Nr+1)';
rbr=mag1*interp1(theev,ybr,theeplotv)+Db/2;
rbmax=max(abs(interp1(theev,ybr,theeplotv)));
if  rbmax~=0
ii=0;
while floor(rbmax)==0
    ii=ii+1;
    rbmax=rbmax*10;
end
end
X=floor(rbmax);
numcircle=2*X+1;
pbmag1=0;
for i=1:length(rbr)
    if rbr(i)<0
        pbmag1=1;
        break;
    end
end
xbore=rbr.*cosd(theeplotv+Tg*180/pi);
ybore=rbr.*sind(theeplotv+Tg*180/pi);
ixbore=Db/2*cosd(theeplotv);
iybore=Db/2*sind(theeplotv);
theeplotv=linspace(dtg,360-dtg,Nplot+1);
theev=linspace(dtg,360-dtg,Nr+1)';
rbr=mag1*interp1(theev,ybr,theeplotv)+Db/2;
rrr=rbr-mag2*interp1(theev,hminv,theeplotv);
pbmag2=0;
for i=1:length(rrr)
    if rrr(i)<0
        pbmag2=1;
        break;
    end
end
xrr=rrr.*cosd(theeplotv+Tg*180/pi);
yrr=rrr.*sind(theeplotv+Tg*180/pi);
fyp=fyring;
flp=flring;

truncatedradialforce=0;
for i=1:length(fyp)
    if abs(fyp(i))>fstop
        truncatedradialforce=1;
        break;
    end
end
truncatedlinerforce=0;
for i=1:length(fyp)
    if abs(flp(i))>fstop
        truncatedlinerforce=1;
        break;
    end
end

thforce=linspace(dtg,360-dtg,numforce);
rrrforce=interp1(theeplotv,rrr,thforce);

fyp=interp1(theeplotv,fyp,thforce);
flp=interp1(theeplotv,flp,thforce);

rmax=max(abs(rrr));

fyp=fyp.*(abs(fyp)<fstop)+fstop*sign(fyp).*(abs(fyp)>=fstop);
flp=flp.*(abs(flp)<fstop)+fstop*sign(flp).*(abs(flp)>=fstop);

fymax=0;
for i=1:length(fyp)
    if (abs(fyp(i))<fstop) && (abs(fyp(i))>fymax)
        fymax=abs(fyp(i)); 
    end
end
flmax=0;
for i=1:length(flp)
    if (abs(flp(i))<fstop) && (abs(flp(i))>flmax)
        flmax=abs(flp(i)); 
    end
end

fyp=rmax/2/fymax*fyp.*(abs(fyp)<fstop)+rmax*sign(fyp).*(abs(fyp)>=fstop);
flp=rmax/2/fymax*flp.*(abs(flp)<fstop)+rmax*sign(flp).*(abs(flp)>=fstop);

numplot=1000;

figure(4)
clf
set(gcf,'Position',[144 64  1302 916])
hold on
plot(ixbore,iybore,'--','Color',[0 0.8 0])
plot(xbore,ybore,'b','linewidth',1.5)
plot(xrr,yrr,'r','linewidth',1.5)

for i=1:numforce
    
    xi=rrrforce(i)*cosd(thforce(i)+Tg*180/pi);
    xf=(rrrforce(i)-fyp(i))*cosd(thforce(i)+Tg*180/pi);
    xplot=linspace(xi,xf,numplot);
    yi=rrrforce(i)*sind(thforce(i)+Tg*180/pi);
    yf=(rrrforce(i)-fyp(i))*sind(thforce(i)+Tg*180/pi);
    yplot=linspace(yi,yf,numplot);
    plot(xplot,yplot,'b')
    
end
if rbmax~=0
theeplotvl=linspace(0,360,Nplot+1);
for j=1:numcircle
   
    diam=Db/2+5*10^(-ii-1)*mag1*j;
    xplot=diam*cosd(theeplotvl);
    yplot=diam*sind(theeplotvl);
    plot(xplot,yplot,'--','Color',[0 0 0])
    x1=(diam+5*10^(-ii-1)*mag1/4)*cosd(75);
    y1=(diam+5*10^(-ii-1)*mag1/4)*sind(75);
    txt1=num2str(5*10^(-ii-1)*j*1e6);
    txt1ext='\mum';
    txt1=strcat(txt1,txt1ext);
    text(x1,y1,txt1)
    
end
end
if rbmax==0
    diam=Db/2;
end
for i=1:12
    xplot=linspace(0,diam*cosd((i-1)*30),1000);
    yplot=linspace(0,diam*sind((i-1)*30),1000);
    plot(xplot,yplot,'--','Color',[0,0,0])
    x1=diam*cosd((i-1)*30)*1.1;
    y1=diam*sind((i-1)*30)*1.1;
    txt1=num2str((i-1)*30);
    text(x1,y1,txt1)
    
end

x1=-diam*1.1;
y1=-diam*1.2^2*0.75;
txt1=num2str(Tg*180/pi);
a='Gap location: ';
txt1=strcat(a,txt1);
a=' deg';
txt1=strcat(txt1,a);
if truncatedradialforce
    a='Truncated force';
    txt1=strcat(txt1,a);
end
if pbmag1
    a=', Bore not plotted correctly because of the magnification coefficient';
    txt1=strcat(txt1,a);
end
if pbmag2
    a=', Ring not plotted correctly because of the magnification coefficient';
    txt1=strcat(txt1,a);
end
text(x1,y1,txt1);
title('Radial force distribution');
axis equal
legend('Not distorted Bore','Magnified distorted Bore','Ring with mignified deformation')

xlim([-diam*1.25 diam*1.25])
ylim([-diam*1.25 diam*1.25])
 set(gca,'FontSize',16);
set(gca,'YTick',[]);
set(gca,'XTick',[]);

figure(5)
clf
set(gcf,'Position',[144 64  1302 916])
hold on
plot(ixbore,iybore,'--','Color',[0 0.8 0])
plot(xbore,ybore,'b','linewidth',1.5)
plot(xrr,yrr,'r','linewidth',1.5)

for i=1:numforce
    
    xi=rrrforce(i)*cosd(thforce(i)+Tg*180/pi);
    xf=(rrrforce(i)-flp(i))*cosd(thforce(i)+Tg*180/pi);
    xplot=linspace(xi,xf,numplot);
    yi=rrrforce(i)*sind(thforce(i)+Tg*180/pi);
    yf=(rrrforce(i)-flp(i))*sind(thforce(i)+Tg*180/pi);
    yplot=linspace(yi,yf,numplot);
    plot(xplot,yplot,'b')
    
end
if rbmax~=0
theeplotvl=linspace(0,360,Nplot+1);
for j=1:numcircle
   
    diam=Db/2+5*10^(-ii-1)*mag1*j;
    xplot=diam*cosd(theeplotvl);
    yplot=diam*sind(theeplotvl);
    plot(xplot,yplot,'--','Color',[0 0 0])
    x1=(diam+5*10^(-ii-1)*mag1/4)*cosd(75);
    y1=(diam+5*10^(-ii-1)*mag1/4)*sind(75);
    txt1=num2str(5*10^(-ii-1)*j*1e6);
    txt1ext='\mum';
    txt1=strcat(txt1,txt1ext);
    text(x1,y1,txt1)
    
end
end
if rbmax==0
    diam=Db/2;
end
for i=1:12
    xplot=linspace(0,diam*cosd((i-1)*30),1000);
    yplot=linspace(0,diam*sind((i-1)*30),1000);
    plot(xplot,yplot,'--','Color',[0,0,0])
    x1=diam*cosd((i-1)*30)*1.1;
    y1=diam*sind((i-1)*30)*1.1;
    txt1=num2str((i-1)*30);
    text(x1,y1,txt1)
    
end

x1=-diam*1.1;
y1=-diam*1.2^2*0.75;
txt1=num2str(Tg*180/pi);
a='Gap location: ';
txt1=strcat(a,txt1);
a=' deg';
txt1=strcat(txt1,a);
if truncatedlinerforce
    a='Truncated force';
    txt1=strcat(txt1,a);
end
if pbmag1
    a=', Bore not plotted correctly because of the magnification coefficient';
    txt1=strcat(txt1,a);
end
if pbmag2
    a=', Ring not plotted correctly because of the magnification coefficient';
    txt1=strcat(txt1,a);
end
text(x1,y1,txt1);
title('Liner ring force distribution');
axis equal
legend('Not distorted Bore','Magnified distorted Bore','Ring with mignified deformation')

xlim([-diam*1.25 diam*1.25])
ylim([-diam*1.25 diam*1.25])
set(gca,'FontSize',16);

set(gca,'YTick',[]);
set(gca,'XTick',[]);

else

theeplotv=linspace(0,360,Nplot+1);
ixbore=Db/2*cosd(theeplotv);
iybore=Db/2*sind(theeplotv);

theeplotv=linspace(dtg,360-dtg,Nplot+1);
theev=linspace(dtg,360-dtg,Nr+1)';

rrmax=max(abs(interp1(theev,hminv,theeplotv)));
if rrmax~=0
ii=0; 
while floor(rrmax)==0
    ii=ii+1;
    rrmax=rrmax*10;
end
X=floor(rrmax);
end

theeplotv=linspace(dtg,360-dtg,Nplot+1);
rrr=Db/2-mag2*interp1(theev,hminv,theeplotv);
pbmag2=0;
for i=1:length(rrr)
    if rrr(i)<0
        pbmag2=1;
        break;
    end
end
xrr=rrr.*cosd(theeplotv+Tg*180/pi);
yrr=rrr.*sind(theeplotv+Tg*180/pi);
fyp=fyring;
flp=flring;

truncatedradialforce=0;
for i=1:length(fyp)
    if abs(fyp(i))>fstop
        truncatedradialforce=1;
        break;
    end
end
truncatedlinerforce=0;
for i=1:length(fyp)
    if abs(flp(i))>fstop
        truncatedlinerforce=1;
        break;
    end
end

thforce=linspace(dtg,360-dtg,numforce);
rrrforce=interp1(theeplotv,rrr,thforce);

fyp=interp1(theeplotv,fyp,thforce);
flp=interp1(theeplotv,flp,thforce);

rmax=max(abs(rrr));

fyp=fyp.*(abs(fyp)<fstop)+fstop*sign(fyp).*(abs(fyp)>=fstop);
flp=flp.*(abs(flp)<fstop)+fstop*sign(flp).*(abs(flp)>=fstop);

fymax=0;
for i=1:length(fyp)
    if (abs(fyp(i))<fstop) && (abs(fyp(i))>fymax)
        fymax=abs(fyp(i)); 
    end
end
flmax=0;
for i=1:length(flp)
    if (abs(flp(i))<fstop) && (abs(flp(i))>flmax)
        flmax=abs(flp(i)); 
    end
end

fyp=rmax/2/fymax*fyp.*(abs(fyp)<fstop)+rmax*sign(fyp).*(abs(fyp)>=fstop);
flp=rmax/2/fymax*flp.*(abs(flp)<fstop)+rmax*sign(flp).*(abs(flp)>=fstop);

numplot=1000;

figure(4)
clf
set(gcf,'Position',[144 64  1302 916])
hold on
plot(ixbore,iybore,'b','linewidth',1.5)
plot(xrr,yrr,'r','linewidth',1.5)

for i=1:numforce
    
    xi=rrrforce(i)*cosd(thforce(i)+Tg*180/pi);
    xf=(rrrforce(i)-fyp(i))*cosd(thforce(i)+Tg*180/pi);
    xplot=linspace(xi,xf,numplot);
    yi=rrrforce(i)*sind(thforce(i)+Tg*180/pi);
    yf=(rrrforce(i)-fyp(i))*sind(thforce(i)+Tg*180/pi);
    yplot=linspace(yi,yf,numplot);
    plot(xplot,yplot,'b')
    
end

if rrmax~=0
numcircle=X+1; 
theeplotvl=linspace(0,360,Nplot+1);
for j=1:numcircle
   
    diam=Db/2-10^(-ii)*mag1*j;
    xplot=diam*cosd(theeplotvl);
    yplot=diam*sind(theeplotvl);
    plot(xplot,yplot,'--','Color',[0 0 0])
    x1=(diam+10^(-ii)*mag1/4)*cosd(75);
    y1=(diam+10^(-ii)*mag1/4)*sind(75);
    txt1=num2str(10^(-ii)*j*1e6);
    txt1ext='\mum';
    txt1=strcat(txt1,txt1ext);
    text(x1,y1,txt1)
    
end
end
diam=Db/2;
for i=1:12
    xplot=linspace(0,diam*cosd((i-1)*30),1000);
    yplot=linspace(0,diam*sind((i-1)*30),1000);
    plot(xplot,yplot,'--','Color',[0,0,0])
    x1=diam*cosd((i-1)*30)*1.1;
    y1=diam*sind((i-1)*30)*1.1;
    txt1=num2str((i-1)*30);
    text(x1,y1,txt1)
    
end

x1=-diam*1.1;
y1=-diam*1.2^2*0.75;
txt1=num2str(Tg*180/pi);
a='Gap location: ';
txt1=strcat(a,txt1);
a=' deg';
txt1=strcat(txt1,a);
if truncatedradialforce
    a='Truncated force';
    txt1=strcat(txt1,a);
end

if pbmag2
    a=', Ring not plotted correctly because of the magnification coefficient';
    txt1=strcat(txt1,a);
end

text(x1,y1,txt1);
title('Radial force distribution');
axis equal
legend('Bore','Ring with mignified deformation')

xlim([-diam*1.15^2 diam*1.15^2])
ylim([-diam*1.15^2 diam*1.15^2])
 set(gca,'FontSize',16);
set(gca,'YTick',[]);
set(gca,'XTick',[]);
figure(5)
clf
set(gcf,'Position',[144 64  1302 916])
hold on
plot(ixbore,iybore,'b','linewidth',1.5)
plot(xrr,yrr,'r','linewidth',1.5)

for i=1:numforce
    
    xi=rrrforce(i)*cosd(thforce(i)+Tg*180/pi);
    xf=(rrrforce(i)-flp(i))*cosd(thforce(i)+Tg*180/pi);
    xplot=linspace(xi,xf,numplot);
    yi=rrrforce(i)*sind(thforce(i)+Tg*180/pi);
    yf=(rrrforce(i)-flp(i))*sind(thforce(i)+Tg*180/pi);
    yplot=linspace(yi,yf,numplot);
    plot(xplot,yplot,'b')
    
end
if rrmax~=0
numcircle=X+1;
theeplotvl=linspace(0,360,Nplot+1);
for j=1:numcircle
   
    diam=Db/2-10^(-ii)*mag1*j;
    xplot=diam*cosd(theeplotvl);
    yplot=diam*sind(theeplotvl);
    plot(xplot,yplot,'--','Color',[0 0 0])
    x1=(diam+10^(-ii)*mag1/4)*cosd(75);
    y1=(diam+10^(-ii)*mag1/4)*sind(75);
    txt1=num2str(10^(-ii)*j*1e6);
    txt1ext='\mum';
    txt1=strcat(txt1,txt1ext);
    text(x1,y1,txt1)
    
end
end
diam=Db/2;
for i=1:12
    xplot=linspace(0,diam*cosd((i-1)*30),1000);
    yplot=linspace(0,diam*sind((i-1)*30),1000);
    plot(xplot,yplot,'--','Color',[0,0,0])
    x1=diam*cosd((i-1)*30)*1.1;
    y1=diam*sind((i-1)*30)*1.1;
    txt1=num2str((i-1)*30);
    text(x1,y1,txt1)
    
end

x1=-diam*1.1;
y1=-diam*1.2^2*0.75;
txt1=num2str(Tg*180/pi);
a='Gap location: ';
txt1=strcat(a,txt1);
a=' deg';
txt1=strcat(txt1,a);
if truncatedlinerforce
    a='Truncated force';
    txt1=strcat(txt1,a);
end

if pbmag2
    a=', Ring not plotted correctly because of the magnification coefficient';
    txt1=strcat(txt1,a);
end
text(x1,y1,txt1);
title('Liner ring force distribution');
axis equal
legend('Bore','Ring with mignified deformation')

xlim([-diam*1.15^2 diam*1.15^2])
ylim([-diam*1.15^2 diam*1.15^2])
 set(gca,'FontSize',16);
set(gca,'YTick',[]);
set(gca,'XTick',[]);
end

end