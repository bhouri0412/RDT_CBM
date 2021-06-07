function [fs,rforce]=rdt_from_ov_to_fs_and_pr(inputs)

%this code computes the radial linear force: 1st column contains angles in
%degrees in ring fixed reference and second one force par circumference 
%length in N/m. it also computes free-shape in usual representation: 
%1st column contains angles in degrees in ring fixed reference and second 
%one the coordintes in millimeters  
%as inputs it needs a text file: 1st column contains angles in degrees in
%ring fixed reference and second one contains the ovality in microns 
%in usual representation

addpath(genpath(pwd))

prompt='Bore diameter (mm): ';
Db=input(prompt)*1e-3;

prompt='Ring Young modulus (GPa): ';
Er=input(prompt)*1e9;

prompt='Gap size when ring is closed to ovality (under constant pressure) (mm): ';
gap=input(prompt)*1e-3;

prompt='Gap size when ring is closed within the bore (mm): ';
gap2=input(prompt)*1e-3;

auo=inputs.ring.auo;
aui=inputs.ring.aui;
alo=inputs.ring.alo;
ali=inputs.ring.ali;
huo=inputs.ring.huo;
hui=inputs.ring.hui;
hlo=inputs.ring.hlo;
hli=inputs.ring.hli;
thrt=inputs.ring.thrt;           
thrb=inputs.ring.thrb;           
rb1=inputs.ring.rb1;
rb2=inputs.ring.rb2;
rbn=inputs.ring.rbn;
a10=inputs.ring.a10;
a11=inputs.ring.a11;
a12=inputs.ring.a12;
a20=inputs.ring.a20;
a21=inputs.ring.a21;
a22=inputs.ring.a22;
arm=inputs.ring.arm;
Ac=inputs.ring.Ac;
Izz=inputs.ring.Izz;
Iyy=inputs.ring.Iyy;
Izr=inputs.ring.Izr;
Ip=inputs.ring.Ip;
Jt=inputs.ring.Jt;
alp=inputs.ring.alp;

Rr=Db/2-arm;
Tgg=gap/2/Rr;

prompt='Number of elements for the FEM: ';
Nbe=input(prompt); 
Nbnod=Nbe+1;
Telt=(2*pi-2*Tgg)/(Nbe);
prompt='Number of points within one element: ';
Npe=input(prompt);
Nr=Nbe*Npe; 

Tr=linspace(Tgg,2*pi-Tgg,Nbe*Npe+1)';

%ask for file containing ovality in usual representation in ring fixed 
%reference: first column containins angles in degrees and second one 
%contains ovality in usual representation in microns
ovdata=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\Theo_OV_U.txt'); %user has to provide this

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

Le=(2*pi*Rr-gap)/Nbe; 
ds=Le/Npe;
Tgg=gap/2/Rr;

kov=((Rr+ov(1:Nr/2+1)).^2+2*dov(1:Nr/2+1).^2-(Rr+ov(1:Nr/2+1)).*d2ov(1:Nr/2+1))./((Rr+ov(1:Nr/2+1)).^2+dov(1:Nr/2+1).^2).^(3/2);

prompt='Constant linear force i.e. pressure times the ring height (axial width) (N/m): ';
Pov=input(prompt);

Tp=linspace(Tgg,pi,Nr/2+1)'; 
Pm=Pov*ones(length(Tp),1);

M=Ring_moment(Tp,Pm,Tr,ov,dov,Tg,Rr);
kapfs=kov-M/Er/Izr;
Kapfs=[kapfs(1:end-1);flipud(kapfs)];

s_span=[pi*Rr:-ds:Tg*Rr]';    

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
Rfs=[flipud(rfs);rfs(2:end)];
fs=[Thefs*180/pi Rfs*1e3];

function [feq]=fseqn(s,y,yp)
    feq=zeros(3,1);
    M_ode=ringM(s,Pm,Tp,Tg,ov,dov,Rr,Tr,Npe,Nbe);
    K_ode=ringK(s,ov,dov,d2ov,Rr,Npe,Nbe);
    feq(1)=yp(1)-sqrt(1-yp(2)^2)/y(2);
    feq(2)=yp(2)-y(3)*yp(1);
    feq(3)=K_ode-M_ode/Er/Izr-(y(2)^2+2*y(3)^2-y(2)*yp(3)*sqrt(y(2)^2+y(3)^2))/(y(2)^2+y(3)^2)^(3/2);

end

%========================bore input==============================%
prompt='Type 1 to take into account bore thermal distortion, 0 otherwise: ';
IsBDist=input(prompt);
if IsBDist
    
    prompt='Type 1 to provide the amplitudes and phases directly and 0 to provide dr distribution: ';
    Ampdata=input(prompt);
    if Ampdata
        prompt='Highest order in the Fourier series for bore distortion: ';
        Norder=input(prompt);
        Ab=zeros(1,Norder);
        if Norder>1
            Phib=zeros(1,Norder-1);
        else
            Phib=0;
        end
        for i=1:Norder
            if i==1
                prompt='Magnitude of 0th order (\mu m): ';
                Ab(1)=input(prompt)*1e-6;
            end
            if i==2
                prompt='Magnitude of 2nd order (\mu m): ';
                Ab(2)=input(prompt)*1e-6;
                prompt='Phase of 2nd order (rad): ';
                Phib(1)=input(prompt);
            end
            if i==3
                prompt='Magnitude of 3rd order (\mu m): '; 
                Ab(3)=input(prompt)*1e-6;
                prompt='Phase of 3rd order (rad): '; 
                Phib(2)=input(prompt);
            end
            if i>3
                strorder=num2str(i);
                s1='Magnitude of';
                s2='th order (\mu m): ';
                prompt=[s1,' ',strorder,s2]; 
                Ab(i)=input(prompt)*1e-6; 
                sp1='Phase of';
                sp2='th order (rad); ';
                prompt=[sp1,' ',storder,sp2]; 
                Phib(i-1)=input(prompt);
            end
        
        end
    else
        
        prompt='Highest order in the Fourier series for bore distortion: ';
        Norder=input(prompt);
        
        %ask for file containing bore distortion: first column contains
        %angles in degrees in bore fixed reference from thrust to thrust 
        %side and second column contains dr in microns
        dr=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\bore_dr.txt'); %user has to provide this
        thdr=dr(:,1);  
        dr=dr(:,2); 
        dr=dr'*1e-6;
        thdr=thdr';
        
        Ab=zeros(1,Norder);
            
        if Norder>1
            Phib=zeros(1,Norder-1);
        else
            Phib=0;
        end
            
        dtheta=2*pi/1000;
        thetaint=linspace(0,2*pi,1001);
        dr=interp1(thdr*pi/180,dr,thetaint,'linear','extrap');
        fint=dr.*cos(pi/180*thetaint);
        ak=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
        fint=dr.*sin(pi/180*thetaint);
        bk=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
        Ab(1)=sqrt(ak^2+bk^2);
            
        if Norder>1
                
           for i=2:Norder
                   
               fint=dr.*cos(i*pi/180*thetaint);
               ak=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
               fint=dr.*sin(i*pi/180*thetaint);
               bk=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
               Ab(i)=sqrt(ak^2+bk^2);
               Phib(i-1)=atan2(ak,bk);
           end
                
        end
            
    end
        
else
    Ab=zeros(1,2); Phib=0;
end

uref=5e-6;  
href=20e-6;

%============================ring input=======================%
prompt='Ring density (kg/m^3): ';
rhor=input(prompt);

prompt='Ring Poisson Ratio: ';
nur=input(prompt);

Gr=Er/(2*(1+nur)); 

prompt='Ring thermal expansion coefficient: ';
alT=input(prompt);

prompt='Ring gap position (deg): ';
Tg=input(prompt)*pi/180;

prompt='Ring temperatrue (C): ';
rTemp=input(prompt); 

Rr=Rr*(1+alT*(rTemp-25));   

Tgg2=gap2/2/Rr;
Telt2=(2*pi-2*Tgg2)/(Nbe); 
Kapfs=interp1(linspace(Tg,2*pi-Tg,Nr+1)',Kapfs,linspace(Tgg2,2*pi-Tgg2,Nr+1)','linear','extrap');
Tgg=Tgg2;
Telt=Telt2;

Tnod=(Tgg:(Telt):2*pi-Tgg)+Tg;

Le=Rr*Telt;
Vnod=[Rr*cos(Tnod),Rr*sin(Tnod)];

mr=rhor*2*pi*Rr*Ac;

%======================liner input================%
prompt='Liner Young modulus (GPa): ';
El=input(prompt)*1e9;

prompt='Liner Poisson Ratio: ';
nul=input(prompt);

prompt='Plateau Ratio: ';
PR=input(prompt); 

prompt='Liner surface roughness standard deviation (\mu m): ';
sigmap=input(prompt)*1e-6; 

%======================piston input================%
prompt='Upper land reference diameter (mm): ';
Drldu=input(prompt)*1e-3; 

prompt='Ring groove root diameter (mm): ';
Drg=input(prompt)*1e-3; 

prompt='Lower land reference diameter (mm): ';
Drldl=input(prompt)*1e-3;

prompt='Ring groove inner axial height (mm): ';
hgi=input(prompt)*1e-3; 

prompt='Ring groove upper flank angle (rad): ';
thgt=input(prompt); 

prompt='Ring groove lower flank angle (rad): ';
thgb=input(prompt);

prompt='Ring groove upper flank surface roughness standard deviation (\mu m): ';
sigmagt=input(prompt)*1e-6;

prompt='Ring groove lower flank surface roughness standard deviation (\mu m): ';
sigmagb=input(prompt)*1e-6; 

prompt='Groove Poisson Ratio: ';
nug=input(prompt); 

prompt='Groove Young modulus (GPa): ';
Eg=input(prompt)*1e9; 

prompt='Radius increase at the upper land during engine operation (\mu m): ';
exp_land_u=input(prompt)*1e-6; 

prompt='Radius increase at the lower land during engine operation (\mu m): ';
exp_land_l=input(prompt)*1e-6; 

prompt='Top ring groove radius increase during engine operation (\mu m): ';
gr_exp=input(prompt)*1e-6;

%======================contact input================%
prompt='z coefficient for the asperity ring/liner and ring/groove contact interaction (6.804 is the adopted value in the simplified formulation): ';
z=input(prompt); 
prompt='K coefficient for the asperity ring/liner and ring/groove contact interaction (1.198e-4 is the adopted value in the simplified formulation): ';
K=input(prompt); 
prompt='A coefficient for the asperity ring/liner and ring/groove contact interaction (4.4068e-5 is the adopted value in the simplified formulation): ';
A=input(prompt);
prompt='Omega coefficient for the asperity ring/liner and ring/groove contact interaction (4 is the adopted value in the simplified formulation): ';
omega=input(prompt); 
prompt='Friction coefficient for the asperity ring/liner and ring/groove contact interaction: ';
fc_dry=input(prompt); 
prompt='Constant coefficient for the asperity ring/liner and ring/groove contact interaction (1 is the adopted value in the simplified formulation): ';
cfct=input(prompt); 
Pk_rl=cfct*K*A*2./((1-nul^2)/El+(1-nul^2)/Er);
Pk_rg=cfct*K*A*2./((1-nug^2)/Eg+(1-nur^2)/Er);

y1u=-aui;
y2u=arm-max((arm-auo),(Db-Drldu)/2+Ab(1)-exp_land_u);
lfu=y2u-y1u;
y1l=-ali;
y2l=arm-max((arm-alo),(Db-Drldl)/2+Ab(1)-exp_land_l);
lfl=y2l-y1l;

%======================groove input================%
prompt='Type 1 to take into account groove thermal distortion and 0 instead: ';
IsGDist=input(prompt);
if IsGDist
    prompt='Type 1 to provide the amplitudes and phases directly for the upper groove deformation or 0 to provide the axial coordinate distribution: ';
    Ampdatau=input(prompt);
    if Ampdatau
        prompt='Number of orders starting from the 2nd one in the Fourier series for the upper groove deformation: ';
        Norderu=input(prompt);
        Agu=zeros(1,Norderu);
        Phigu=zeros(1,Norderu);
        for i=1:Norderu
            strorder=num2str(i+1);
            s1='Magnitude of';
            s2='th order (\mu m): ';
            prompt=[s1,' ',strorder,s2]; 
            Agu(i)=input(prompt)*1e-6; 
            sp1='Phase of';
            sp2='th order (rad); ';
            prompt=[sp1,' ',storder,sp2];
            Phigu(i)=input(prompt);
        end
    else
        prompt='Number of orders starting from the 2nd one in the Fourier series for the upper groove deformation: ';
        Norderu=input(prompt);
        
        %ask for file containing groove upper flank distortion: first 
        %column contains angles in degrees in bore fixed reference from 
        %thrust to thrust side and second column contains dz in microns
        dzgu=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\upper_flank_dz.txt'); %user has to provide this
        thdrgu=dr(:,1);  
        dzgu=dzgu(:,2); 
        
        dzgu=dzgu'*1e-6;
        thdrgu=thdrgu';
        
        Agu=zeros(1,Norderu);
        Phigu=zeros(1,Norderu);
        
        dtheta=2*pi/1000;
        thetaint=linspace(0,2*pi,1001);
        dzgu=interp1(thdrgu*pi/180,dzgu,thetaint,'linear','extrap');
        
        for i=1:Norderu
            
            fint=dzgu.*cos((i+1)*pi/180*thetaint);
            ak=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
            fint=dzgu.*sin((i+1)*pi/180*thetaint);
            bk=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
            Agu(i)=sqrt(ak^2+bk^2);
            Phigu(i)=atan2(ak,bk);
        
        end
        
    end
    prompt='Type 1 to provide the amplitudes and phases directly for the lower groove deformation or 0 to provide the axial coordinate distribution: ';
    Ampdatal=input(prompt);
    if Ampdatal
        prompt='Number of orders starting from the 2nd one in the Fourier series for the lower groove deformation: ';
        Norderl=input(prompt);
        Agl=zeros(1,Norderl);
        Phigl=zeros(1,Norderl);
        for i=1:Norderl
            strorder=num2str(i+1);
            s1='Magnitude of';
            s2='th order (\mu m): ';
            prompt=[s1,' ',strorder,s2]; 
            Agl(i)=input(prompt)*1e-6; 
            sp1='Phase of';
            sp2='th order (rad); ';
            prompt=[sp1,' ',storder,sp2];
            Phigl(i)=input(prompt);
        end
    else
        prompt='Number of orders starting from the 2nd one in the Fourier series for the lower groove deformation: ';
        Norderl=input(prompt);
        
        %ask for file containing groove lower flank distortion: first 
        %column contains angles in degrees in bore fixed reference from 
        %thrust to thrust side and second column contains dz in microns
        dzgl=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\lower_flank_dz.txt'); %user has to provide this
        thdrgl=dr(:,1);  
        dzgl=dzgl(:,2); 
        
        dzgl=dzgl'*1e-6;
        thdrgl=thdrgl';
        
        Agl=zeros(1,Norderl);
        Phigl=zeros(1,Norderl);
        
        dtheta=2*pi/1000;
        thetaint=linspace(0,2*pi,1001);
        dzgl=interp1(thdrgl*pi/180,dzgl,thetaint,'linear','extrap');
        
        for i=1:Norderl
            
            fint=dzgl.*cos((i+1)*pi/180*thetaint);
            ak=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
            fint=dzgl.*sin((i+1)*pi/180*thetaint);
            bk=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
            Agl(i)=sqrt(ak^2+bk^2);
            Phigl(i)=atan2(ak,bk);
        
        end
        
    end
   
else
    Agu=[0,0,0]*0;
    Phigu=[0,0,0]; 
    Agl=[0,0,0]*0; 
    Phigl=[0,0,0];
end

prompt='Groove thermal tilting - upper flank (rad): '; 
tilt_thu=input(prompt);

prompt='Groove thermal tilting - lower flank (rad): '; 
tilt_thl=input(prompt);

hrc=hui+thrt*aui+hli+thrb*ali;
hgc=hgi+(Db/2-Drg/2-gr_exp-arm)*(thgt+thgb);
gcl=hgc-hrc;

%=====================prepare local input matrix======================%
inputsl2.bore.Db=Db;
inputsl2.bore.Ab=Ab;
inputsl2.bore.Phib=Phib;
inputsl2.ring.auo=auo;
inputsl2.ring.aui=aui;
inputsl2.ring.alo=alo;
inputsl2.ring.ali=ali;
inputsl2.ring.huo=huo;
inputsl2.ring.hui=hui;
inputsl2.ring.hlo=hlo;
inputsl2.ring.hli=hli;
inputsl2.ring.thrt=thrt;           
inputsl2.ring.thrb=thrb;           
inputsl2.ring.rb1=rb1;
inputsl2.ring.rb2=rb2;
inputsl2.ring.rbn=rbn;
inputsl2.ring.a10=a10;
inputsl2.ring.a11=a11;
inputsl2.ring.a12=a12;
inputsl2.ring.a20=a20;
inputsl2.ring.a21=a21;
inputsl2.ring.a22=a22;
inputsl2.ring.arm=arm;
inputsl2.ring.gap=gap;
inputsl2.ring.rhor=rhor;
inputsl2.ring.Er=Er;
inputsl2.ring.alT=alT;
inputsl2.ring.nur=nur;
inputsl2.ring.Gr=Gr;
inputsl2.ring.Ac=Ac;
inputsl2.ring.Izz=Izz;
inputsl2.ring.Iyy=Iyy;
inputsl2.ring.Izr=Izr;
inputsl2.ring.Ip=Ip;
inputsl2.ring.Jt=Jt;
inputsl2.ring.alp=alp;
inputsl2.ring.Rr=Rr;
inputsl2.ring.Le=Le;
inputsl2.ring.mr=mr;
inputsl2.ring.Tg=Tg;
inputsl2.ring.rTemp=rTemp;

inputsl2.piston.Drldu=Drldu;
inputsl2.piston.Drg=Drg;
inputsl2.piston.Drldl=Drldl;
inputsl2.piston.hgi=hgi;
inputsl2.piston.thgt=thgt;
inputsl2.piston.thgb=thgb;
inputsl2.piston.sigmag1t=sigmagt;
inputsl2.piston.sigmag1b=sigmagb;
inputsl2.piston.nug=nug;
inputsl2.piston.Eg=Eg;
inputsl2.piston.exp_land_u=exp_land_u;
inputsl2.piston.exp_land_l=exp_land_l;
inputsl2.piston.gr_exp=gr_exp;

inputsl2.liner.El=El;
inputsl2.liner.nul=nul;
inputsl2.liner.sigmap=sigmap;
inputsl2.liner.PR=PR;

inputsl2.FEM.Nbe=Nbe;
inputsl2.FEM.Nbnod=Nbnod;
inputsl2.FEM.Telt=Telt;
inputsl2.FEM.Npe=Npe;
inputsl2.FEM.Nr=Nr;
inputsl2.FEM.Le=Le;
inputsl2.FEM.Tnod=Tnod;
inputsl2.FEM.Vnod=Vnod;
inputsl2.FEM.uref=uref;
inputsl2.FEM.href=href;

inputsl2.contact.z=z;
inputsl2.contact.Pk_rl=Pk_rl;
inputsl2.contact.Pk_rg=Pk_rg;
inputsl2.contact.omega=omega;
inputsl2.contact.fc_dry=fc_dry;
inputsl2.contact.cfct=cfct;

inputsl2.groove.y1u=y1u;
inputsl2.groove.y2u=y2u;
inputsl2.groove.lfu=lfu;
inputsl2.groove.y1l=y1l;
inputsl2.groove.y2l=y2l;
inputsl2.groove.lfl=lfl;
inputsl2.groove.gcl=gcl;
inputsl2.groove.Agu=Agu;
inputsl2.groove.Phigu=Phigu;
inputsl2.groove.Agl=Agl;
inputsl2.groove.Phigl=Phigl;
inputsl2.groove.tilt_thu=tilt_thu;
inputsl2.groove.tilt_thl=tilt_thl;

%=====================solver======================%
Ug1=zeros(8*Nbnod,1);
Ug1(7:8:end-1)=alp;

F1=Finitial_fs_o(Ug1,Er,Gr,Iyy,Izz,Jt,alp,Rr,Telt,Nbe,Npe,Kapfs);

inputsl2.ring.F1=F1;

[Kg_nd]=Kg_ring_curved_nd_o(Er,Gr,Iyy,Jt,Rr,Telt,Nbe);

inputsl2.FEM.Kg_nd=Kg_nd;

UgIC_nd=zeros(8*Nbnod,1);
[ybn]=bore_dist_static_o(inputsl2);
Tnodd=(0:2*pi/Nbe:2*pi)'+Tg;

id1=1:3:3*Nbnod-2;
id2=2:3:3*Nbnod-1;
id3=3:3:3*Nbnod;

ybn1=ybn(id1);
ybn2=ybn(id2);
ybn3=ybn(id3);

ybn(id1)=interp1(Tnodd,ybn1,Tnod');
ybn(id2)=interp1(Tnodd,ybn2,Tnod');
ybn(id3)=interp1(Tnodd,ybn3,Tnod');

hl_ig=(omega-((abs(F1(1:8:end-7))*2/(rb1+rb2)/Le)/Pk_rl/PR).^(1/z))*sigmap;
UgIC_nd(1:8:8*Nbnod-7)=ybn(id1)-hl_ig/uref;
UgIC_nd(2:8:8*Nbnod-6)=ybn(id2);
UgIC_nd(3:8:8*Nbnod-5)=ybn(id3);

dz_g=(Rr*cos(Tnod')-0)*0;
[zgnu,~]=groove_dist_static_o(inputsl2);
zgnu1=zgnu(id1);
zgnu2=zgnu(id2);
zgnu3=zgnu(id3);

zgnu(id1)=interp1(Tnodd,zgnu1,Tnod');
zgnu(id2)=interp1(Tnodd,zgnu2,Tnod');
zgnu(id3)=interp1(Tnodd,zgnu3,Tnod');

hg_ig=(omega-((abs(F1(4:8:end-4))*2/Le/lfu)/Pk_rg).^(1/z))*sigmagt;
UgIC_nd(4:8:8*Nbnod-4)=(gcl/2-hg_ig+dz_g)/href+zgnu(id1);
UgIC_nd(5:8:8*Nbnod-3)=-Rr*sin(Tnod')*0*Telt/href+zgnu(id2);
UgIC_nd(6:8:8*Nbnod-2)=-Rr*cos(Tnod')*0*Telt^2/href+zgnu(id3);

beta=0*cos(Tnod')+tilt_thu; 
UgIC_nd(7:8:8*Nbnod-1)=(beta)/href*Rr;
UgIC_nd(8:8:8*Nbnod)=-0*sin(Tnod')/href*Rr*Telt;

Ug_nd=UgIC_nd;

prompt='Tolerance for Newton-Raphson algorithm convergence (suggested value: 1e-6): ';
Newton_tol=input(prompt); 
prompt='Maximum number of Newton-Raphson algorithm iterations (suggested value: 100): ';
ItMax=input(prompt);
judge1=1;
judge2=1;
NLS_0=Static_Error_o(Ug_nd,F1,inputsl2);
NLS_k=NLS_0;
j=0;

while j<=ItMax&&(judge1||judge2)
    j=j+1;
    [Ja_stac]=Jacobian_Estac_o(Ug_nd,inputsl2);
    Nstep=-Ja_stac\NLS_k;
    Ugo_nd=Ug_nd;
    Erro=norm(NLS_k);
    itGB=0; dErr=-1;
    
    if j<2
        while (dErr<0)&&(itGB<20)
            if itGB>0
                Nstep=Nstep/4;
            end
            Ug_nd=Ugo_nd+Nstep; 
            NLS_k=Static_Error_o(Ug_nd,F1,inputsl2);
            Err=norm(NLS_k);
            dErr=Erro-Err;
            itGB=itGB+1;
        end
    else
        while (dErr<0)&&(itGB<10)
            if itGB>0
                Nstep=Nstep/2;
            end
            Ug_nd=Ugo_nd+Nstep; 
            NLS_k=Static_Error_o(Ug_nd,F1,inputsl2);
            Err=norm(NLS_k);
            dErr=Erro-Err;
            itGB=itGB+1;
        end
    end

    judge1=norm(Nstep)>Newton_tol;
    judge2=norm(NLS_k)>Newton_tol;
    
    norm_NLS_k(j)=norm(NLS_k);
    norm_step(j)=norm(Nstep);
    
    figure(10)
    cla
    subplot(211)
    semilogy(norm_NLS_k)
    grid
    subplot(212)
    semilogy(norm_step)
    grid
    pause(0.0001)
    Err_NLS_k=norm(NLS_k)
    Err_Nstep=norm(Nstep)

end

dx=1/Npe;
x=0:dx:1;

N1=1-10*x.^3+15*x.^4-6*x.^5;
N2=x-6*x.^3+8*x.^4-3*x.^5;
N3=(x.^2)/2-3*(x.^3)/2+3*(x.^4)/2-(x.^5)/2;
N4=10*x.^3-15*x.^4+6*x.^5;
N5=-4*x.^3+7*x.^4-3*x.^5;
N6=(x.^3)/2-x.^4+x.^5/2;
N=[N1;N2;N3;N4;N5;N6];

Nt1=1-3*x.^2+2*x.^3;
Nt2=x-2*x.^2+x.^3;
Nt3=3*x.^2-2*x.^3;
Nt4=-x.^2+x.^3;

Nt=[Nt1;Nt2;Nt3;Nt4];

[ybn]=bore_dist_static_o(inputsl2);
Tnodd=(0:2*pi/Nbe:2*pi)'+Tg;
id1=1:3:3*Nbnod-2;
id2=2:3:3*Nbnod-1;
id3=3:3:3*Nbnod;

ybn1=ybn(id1);
ybn2=ybn(id2);
ybn3=ybn(id3);

ybn(id1)=interp1(Tnodd,ybn1,Tnod');
ybn(id2)=interp1(Tnodd,ybn2,Tnod');
ybn(id3)=interp1(Tnodd,ybn3,Tnod');

fclr=zeros(1,Nr+1);
ybr=zeros(1,Nr+1);
yrr=zeros(1,Nr+1);
for j=1:Nbe
    ind=1+(j-1)*Npe:1+j*Npe;
    
    Uey=Ug_nd([8*(j-1)+1:8*(j-1)+3,8*j+1:8*j+3]);
    Uet=Ug_nd([8*j-1:8*j,8*j+7:8*j+8]);
    
    yr=(Uey')*N*uref;
    alr=(Uet')*Nt*href/Rr;

    id=3*(j-1)+1:3*(j+1);
    yb=(ybn(id))'*N*uref;
    
    ybr(ind)=yb;
    yrr(ind)=yr;
    
    Nbz=100;
    dz=(rb1+rb2)/Nbz;
    zz=linspace(-rb1,rb2,Nbz+1);
    zz=repmat(zz',1,Npe+1)+rbn;
    
    delybr=repmat(yb-yr,Nbz+1,1);
    alr2=repmat(alr,Nbz+1,1);
    
    ringface=(a12*(zz-rbn).^2+a11*(zz-rbn)+a10).*(zz<=rbn)+(a22*(zz-rbn).^2+a21*(zz-rbn)+a20).*(zz>rbn);
    
    hh=delybr+ringface+alr2.*zz;
    Pcl=PR*Pk_rl*(omega-hh/sigmap).^z.*(hh/sigmap<=omega);
    fcl=sum((Pcl(1:end-1,:)+Pcl(2:end,:))/2*dz);
        
    fclr(ind)=fcl;

end

thee=linspace(Tgg*180/pi,360-Tgg*180/pi,Nr+1);
display(['Ft(CBM) = ', num2str(sum(fclr*Rr/Nr)),'N']) %tengential load

rforce=[thee' fclr'];

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

figure(2)
clf;
hold on
plot(thee,fclr,'LineWidth',1.5)
plot((Tnod-Tg)/pi*180,fclr(1:Npe:end),'ro')
xlabel('Circumferential Direction (degree)')
ylabel('Ring Liner Force (N/m)','FontSize',12)
set(gca,'XTick',[0:90:360])
hold off
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