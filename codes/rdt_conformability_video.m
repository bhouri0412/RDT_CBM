function [res_Tg,res_theer,res_hmin,res_zr,res_yr,res_alr,res_z0,res_uoc,res_uic,res_loc,res_lic,res_ybr,...
    res_thee,res_fy,res_fl,res_fz,res_m,res_Mfinal,res_S_uo,res_S_ui,res_S_lo,res_S_li,util_data]=rdt_conformability_video(inputs)

%this code performs the calculation for the conformability study for different ring gap positions equally spaced between 0 and 359 degrees

addpath(genpath(pwd))

%=============features================================%
prompt='Type 1 to take into account bore thermal distortion, 0 otherwise: ';
IsBDist=input(prompt);

prompt='Type 1 to take into account piston dynamic tilt, 0 otherwise: ';
IsPTilt=input(prompt);

prompt='Type 1 to take into account groove thermal distortion, 0 otherwise: ';
IsGDef=input(prompt);

prompt='Type 1 to provide free shape or its curvature and 0 to provide radial pressure distribution for the circular shape: ';
IsFreeshape=input(prompt);

inputsl4.features.IsBDist=IsBDist;
inputsl4.features.IsPTilt=IsPTilt;
inputsl4.features.IsGDef=IsGDef;
inputsl4.features.IsFreeshape=IsFreeshape;

%============================ring input=======================%
prompt='Bore diameter (mm): ';
Db=input(prompt)*1e-3;
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
prompt='Gap size when ring is closed within the bore (mm): ';
gap=input(prompt)*1e-3; 
Ac=inputs.ring.Ac; 
Izz=inputs.ring.Izz;
Iyy=inputs.ring.Iyy;
Iyr=inputs.ring.Iyr;
Izr=inputs.ring.Izr;
Ip=inputs.ring.Ip;
Jt=inputs.ring.Jt; 
alp=inputs.ring.alp;

prompt='Ring density (kg/m^3): ';
rhor=input(prompt);

prompt='Ring Young modulus (GPa): ';
Er=input(prompt)*1e9;

prompt='Ring Poisson Ratio: ';
nur=input(prompt); 

Gr=Er/(2*(1+nur)); 

prompt='Ring thermal expansion coefficient (1/K): ';
alT=input(prompt); 

prompt='Ring temperatrue (C): ';
rTemp=input(prompt); 

Rro=Db/2-arm;
Rr=Rro*(1+alT*(rTemp-25));  
Tgg=gap/2/Rr;

%ask for file containing ring temperature increase from id to od :first 
%column contains angles in degrees in ring fixed reference and second 
%column contains radial ring temperature increase from id to od in Celsius
dTemp=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\dTemp.txt'); %user has to provide this
angle_dTemp=dTemp(:,1);
dTemp=dTemp(:,2);

%=========================FEM input============================%
prompt='Number of elements for the FEM: ';
Nbe=input(prompt); 
Nbnod=Nbe+1; 
Telt=(2*pi-2*Tgg)/(Nbe); 
prompt='Number of points within one element for contact and lubrication force calculation: ';
Npe=input(prompt); 
Nr=Nbe*Npe; 
prompt='Number of points within one element for gas force calculation: ';
Nge=input(prompt); 
Ng=Nbe*Nge;

Le=Rr*Telt;
id=(1:Ng+1)';
mr=rhor*2*pi*Rro*Ac;

%========================bore input==============================%
if IsBDist 
    prompt='Type 1 to provide the amplitudes and phases directly and 0 to provide dr distribution: ';
    Ampdata=input(prompt);
    if Ampdata
        prompt='Highest order in the Fourier series for bore distrtion: ';
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

%======================== Pistong Secondary Motion ====================%
if IsPTilt 
    prompt='Piston tilt angle (deg): ';
    beta_p=input(prompt)*pi/180; 
else 
    beta_p=0;
end
prompt='Piston offset (mm): ';
off=input(prompt)*1e-3; 

uref=5e-6; 
href=beta_p*Db/2; 

if href==0
    href=20e-6;
end

%======================== Oil Data ========================%
prompt='zk coefficient for oil viscosity calculation: ';
zk=input(prompt); 
prompt='temp1 coefficient for oil viscosity calculation (C): ';
temp1=input(prompt); 
prompt='temp2 coefficient for oil viscosity calculation (C): ';
temp2=input(prompt); 
prompt='Oil density (kg/m^3): ';
rho_oil=input(prompt); 
prompt='hlratio coefficient for oil viscosity calculation: ';
hlratio=input(prompt);
prompt='Beta1 coefficient for oil viscosity calculation: ';
bta1=input(prompt);
prompt='Beta2 coefficient for oil viscosity calculation: ';
bta2=input(prompt); 
zm0=1;

angle_thetag=linspace(Tgg*180/pi,360-Tgg*180/pi,Ng+1);
angle_theta=linspace(Tgg*180/pi,360-Tgg*180/pi,Nr+1);
anglesg=linspace(0,360,Ng+1);
anglesr=linspace(0,360,Nr+1);

angle_ringg=linspace(Tgg*180/pi,360-Tgg*180/pi,1000*Nbe+1);

dTemp=interp1(angle_dTemp,dTemp,angle_ringg,'linear','extrap');

%ask for file containing oil viscosity in groove: first column contains 
%angles in degrees in bore fixed reference from thrust to thrust side and 
%second column contains oil viscority in the groove in Pa.s
mu_oilglob=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\mu_oil.txt'); %user has to provide this
angle_mu_oil=mu_oilglob(:,1);
mu_oilglob=mu_oilglob(:,2);

mu_oilglob=interp1(angle_mu_oil,mu_oilglob,anglesg,'linear','extrap');

%======================liner input================%

prompt='Liner Young modulus (GPa): ';
El=input(prompt)*1e9; 

prompt='Liner Poisson Ratio: ';
nul=input(prompt); 

prompt='Plateau Ratio: ';
PR=input(prompt); 

prompt='Liner roughness (\mu m): ';
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

prompt='Ring groove upper flank roughness (\mu m): ';
sigmagt=input(prompt)*1e-6; 

prompt='Ring groove lower flank roughness (\mu m): ';
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

if IsGDef
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
        
        dtheta=2*pi/10000;
        thetaint=linspace(0,2*pi,10001);
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
        
        dtheta=2*pi/10000;
        thetaint=linspace(0,2*pi,10001);
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
    Agu=[0,0,0]*1e-6*0;  
    Phigu=[0,0,0]; 
    Agl=[0,0,0]*1e-6*0; 
    Phigl=[0,0,0];
end

prompt='Groove thermal tilting - upper flank (rad): '; 
tilt_thu=input(prompt);

prompt='Groove thermal tilting - lower flank (rad): ';
tilt_thl=input(prompt);

hrc=hui+thrt*aui+hli+thrb*ali;
hgc=hgi+(Db/2-Drg/2-gr_exp-arm)*(thgt+thgb);
gcl=hgc-hrc;

%=====================lubrication&dynamic input======================%
prompt='Minimum oil film thickness to consider for liner lubrication in microns (suggested value: half the liner surface roughness standard deviation): ';
minoil=input(prompt)*1e-6;
prompt='Ratio between oil film thickness on liner by liner surface roughness standard deviation to consider as threshold between parially and fully flooded cases (suggested value: 10): ';
oilthreshold=input(prompt);
prompt='Factor by which the oil viscosity should be divided in regions where we have fuel spots (suggested value: 100): ';
viscosityfactor=input(prompt);

%ask for file containing gas pressure in groove upper region: first column
%contains angles in degrees in bore fixed reference from thrust to thrust 
%side and second column contains gas pressure in groove upper region in Bar
Puglob=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\Pu.txt'); %user has to provide this
angle_Pu=Puglob(:,1);
Puglob=Puglob(:,2);
%ask for file containing gas pressure in groove inner region: first column
%contains angles in degrees in bore fixed reference from thrust to thrust 
%side and second column contains gas pressure in groove inner region in Bar
Piglob=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\Pi.txt'); %user has to provide this
angle_Pi=Piglob(:,1);
Piglob=Piglob(:,2);
%ask for file containing gas pressure in groove lower region: first column
%contains angles in degrees in bore fixed reference from thrust to thrust 
%side and second column contains gas pressure in groove lower region in Bar
Pdglob=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\Pd.txt'); %user has to provide this
angle_Pd=Pdglob(:,1);
Pdglob=Pdglob(:,2);

Puglob=1e5*Puglob;
Piglob=1e5*Piglob;
Pdglob=1e5*Pdglob;

Piglob=interp1(angle_Pi,Piglob,anglesg,'linear','extrap');
Pdglob=interp1(angle_Pd,Pdglob,anglesg,'linear','extrap');
Puglob=interp1(angle_Pu,Puglob,anglesg,'linear','extrap');

prompt='piston speed (m/s): ';
Vp=input(prompt); 
prompt='piston acceleration (m^2/s): ';
Ap=input(prompt);
Fi=-mr*Ap;

%ask for file containing oil film thickness left on liner into the ring: 
%first column contains angles in degrees in bore fixed reference from 
%thrust to thrust side and second column contains oil film thickness left 
%on liner into the ring in microns
hrlsglob=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\hrls.txt'); %user has to provide this
angle_hrls=hrlsglob(:,1);
hrlsglob=hrlsglob(:,2);

hrlsglob=1e-6*hrlsglob;
hrlsglob=interp1(angle_hrls,hrlsglob,anglesr,'linear','extrap');

%ask for file containing oil film thickness on groove upper flank: first 
%column contains angles in degrees in bore fixed reference from thrust to 
%thrust side and second column contains oil film thickness on groove upper
%flank in microns
hotglob=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\hot.txt'); %user has to provide this
angle_hot=hotglob(:,1);
hotglob=hotglob(:,2);

hotglob=1e-6*hotglob;
hotglob=interp1(angle_hot,hotglob,anglesg,'linear','extrap');

%ask for file containing oil film thickness on groove lower flank: first 
%column contains angles in degrees in bore fixed reference from thrust to 
%thrust side and second column contains oil film thickness on groove lower
%flank in microns
hobglob=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\hob.txt'); %user has to provide this
angle_hob=hobglob(:,1);
hobglob=hobglob(:,2);

hobglob=1e-6*hobglob;
hobglob=interp1(angle_hob,hobglob,anglesg,'linear','extrap');

%ask for file containing liner temperature that will be used to compute oil
%viscosity: first column contains angles in degrees in bore fixed reference
%from thrust to thrust side and second column contains liner temperature in
%Celsius
templglob=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\templ.txt'); %user has to provide this
angle_templ=templglob(:,1);
templglob=templglob(:,2);

templglob=interp1(angle_templ,templglob,anglesr,'linear','extrap');

%ask for file containing fuel spots on liner: first column contains angles 
%in degrees in bore fixed reference from thrust to thrust side and second 
%column contains fuel spots (1 for existence of fuel)
isfuelglob=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\isfuel.txt'); %user has to provide this
angle_fuel=isfuelglob(:,1);
isfuelglob=isfuelglob(:,2);

isfuelglob=interp1(angle_fuel,isfuelglob,anglesr,'linear','extrap');

%=====================fully flooded correlation======================%
prompt='b coefficient for fully flooded correlation (2.1): ';
b=input(prompt); 
prompt='c coefficient for fully flooded correlation (0.102): ';
c=input(prompt);
prompt='d coefficient for fully flooded correlation (-0.04): ';
d=input(prompt); 
prompt='e coefficient for fully flooded correlation (0.90229): ';
e=input(prompt); 

%=====================deterministic correlation======================%
prompt='Correlation pressure Pocr for deterministic partially flooded forces (bar) (398.56): ';
Ph=input(prompt)*1e5; 
prompt='Correlation exponent Kocr for deterministic partially flooded forces (2.523): ';
KOCR=input(prompt); 
prompt='Reference muU for deterministic partially flooded forces (Pa.m) (1.5e-2): ';
muUot=input(prompt); 

prompt='ap coefficient for deterministic partially flooded forces (0.039342): ';
ap=input(prompt); 
prompt='cp1 coefficient for exponent Kp in deterministic partially flooded forces (-2.965331): ';
cp1=input(prompt); 
prompt='cp2 coefficient for exponent Kp in deterministic partially flooded forces (1.499148): ';
cp2=input(prompt);
prompt='F0 coefficient for deterministic partially flooded shear force (0.098577): ';
F0=input(prompt); 
prompt='cf1 coefficient for exponent Kp in deterministic partially flooded forces (-0.383954): ';
cf1=input(prompt); 
prompt='cf2 coefficient for exponent Kp in deterministic partially flooded forces (0.138443): ';
cf2=input(prompt); 
prompt='Ring width (mm): ';
rw=input(prompt)*1e-3;

%=====================prepare local input matrix======================%
inputsl4.pistonSM.beta_p=beta_p;
inputsl4.pistonSM.off=off;

inputsl4.bore.Db=Db;
inputsl4.bore.Ab=Ab;
inputsl4.bore.Phib=Phib;

inputsl4.oil.zk=zk;
inputsl4.oil.temp1=temp1;
inputsl4.oil.temp2=temp2;
inputsl4.oil.rho_oil=rho_oil;
inputsl4.oil.hlratio=hlratio;
inputsl4.oil.bta1=bta1;
inputsl4.oil.bta2=bta2;
inputsl4.oil.zm0=zm0;

inputsl4.ring.auo=auo;
inputsl4.ring.aui=aui;
inputsl4.ring.alo=alo;
inputsl4.ring.ali=ali;
inputsl4.ring.huo=huo;
inputsl4.ring.hui=hui;
inputsl4.ring.hlo=hlo;
inputsl4.ring.hli=hli;
inputsl4.ring.thrt=thrt;           
inputsl4.ring.thrb=thrb;           
inputsl4.ring.rb1=rb1;
inputsl4.ring.rb2=rb2;
inputsl4.ring.rbn=rbn;
inputsl4.ring.a10=a10;
inputsl4.ring.a11=a11;
inputsl4.ring.a12=a12;
inputsl4.ring.a20=a20;
inputsl4.ring.a21=a21;
inputsl4.ring.a22=a22;
inputsl4.ring.arm=arm;
inputsl4.ring.gap=gap;
inputsl4.ring.rhor=rhor;
inputsl4.ring.Er=Er;
inputsl4.ring.alT=alT;
inputsl4.ring.nur=nur;
inputsl4.ring.Gr=Gr;
inputsl4.ring.Ac=Ac;
inputsl4.ring.Izz=Izz;
inputsl4.ring.Iyy=Iyy;
inputsl4.ring.Izr=Izr;
inputsl4.ring.Ip=Ip;
inputsl4.ring.Jt=Jt;
inputsl4.ring.alp=alp;
inputsl4.ring.Rr=Rr;
inputsl4.ring.Le=Le;
inputsl4.ring.mr=mr;

inputsl4.ring.rTemp=rTemp;
inputsl4.ring.dTemp=dTemp;

inputsl4.piston.Drldu=Drldu;
inputsl4.piston.Drg=Drg;
inputsl4.piston.Drldl=Drldl;
inputsl4.piston.hgi=hgi;
inputsl4.piston.thgt=thgt;
inputsl4.piston.thgb=thgb;
inputsl4.piston.sigmag1t=sigmagt;
inputsl4.piston.sigmag1b=sigmagb;
inputsl4.piston.nug=nug;
inputsl4.piston.Eg=Eg;
inputsl4.piston.exp_land_u=exp_land_u;
inputsl4.piston.exp_land_l=exp_land_l;
inputsl4.piston.gr_exp=gr_exp;

inputsl4.liner.El=El;
inputsl4.liner.nul=nul;
inputsl4.liner.sigmap=sigmap;
inputsl4.liner.PR=PR;

inputsl4.FEM.Nbe=Nbe;
inputsl4.FEM.Nbnod=Nbnod;
inputsl4.FEM.Telt=Telt;
inputsl4.FEM.Npe=Npe;
inputsl4.FEM.Nr=Nr;
inputsl4.FEM.Le=Le;

inputsl4.FEM.uref=uref;
inputsl4.FEM.href=href;
inputsl4.FEM.Nge=Nge;
inputsl4.FEM.Ng=Ng;

inputsl4.contact.z=z;
inputsl4.contact.Pk_rl=Pk_rl;
inputsl4.contact.Pk_rg=Pk_rg;
inputsl4.contact.omega=omega;
inputsl4.contact.fc_dry=fc_dry;
inputsl4.contact.cfct=cfct;

inputsl4.groove.y1u=y1u;
inputsl4.groove.y2u=y2u;
inputsl4.groove.lfu=lfu;
inputsl4.groove.y1l=y1l;
inputsl4.groove.y2l=y2l;
inputsl4.groove.lfl=lfl;
inputsl4.groove.gcl=gcl;
inputsl4.groove.Agu=Agu;
inputsl4.groove.Phigu=Phigu;
inputsl4.groove.Agl=Agl;
inputsl4.groove.Phigl=Phigl;
inputsl4.groove.tilt_thu=tilt_thu;
inputsl4.groove.tilt_thl=tilt_thl;

inputsl4.engine.Vp=Vp;
inputsl4.engine.Ap=Ap;
inputsl4.ring.Fi=Fi;

inputsl4.liner.PR=PR;
inputsl4.liner.Ph=Ph;
inputsl4.liner.KOCR=KOCR;
inputsl4.liner.ap=ap;
inputsl4.liner.cp1=cp1;
inputsl4.liner.cp2=cp2;
inputsl4.liner.F0=F0;
inputsl4.liner.cf1=cf1;
inputsl4.liner.cf2=cf2;
inputsl4.liner.rw=rw;
inputsl4.liner.muUot=muUot;
inputsl4.liner.b=b;
inputsl4.liner.c=c;
inputsl4.liner.d=d;
inputsl4.liner.e=e;

inputsl4.lubrication.minoil=minoil;
inputsl4.lubrication.oilthreshold=oilthreshold;
inputsl4.lubrication.viscosityfactor=viscosityfactor;

%=============freeshape================%

if IsFreeshape
    prompt='Type 1 to provide the free shape coordinates and 0 to provide the free shape curvature: ';
    IsFreeshaped=input(prompt);
    if IsFreeshaped
        %ask for file containing free-shape coordinates in ring fixed
        %reference: first column containins angles in degrees and second 
        %one contains the freeshape coordinates in usual representation in
        %millimeters
        FS_raw=dlmread('ringfreeshape.txt'); %user has to provide this
        Tri=FS_raw(:,1)*pi/180;
        rfsim=FS_raw(:,2)*1e-3;
        
        Inod=1:1000:1000*Nbe+1;
        dTrfs=(Tri(end)-Tri(1))/(1000*Nbe);
        Trfs=(Tri(1):dTrfs:Tri(length(Tri)))';
        Teltfs=(Trfs(end)-Trfs(1))/Nbe;
        Tnodfs=Trfs(Inod);
    
        rfsi=interp1(Tri,rfsim,Trfs);
        [Ccs]=Cs_3rd_cont_vg(Tnodfs);
        [rnodfs]=ls_rfsi_vg(rfsi,Nbnod,Teltfs,Ccs);
        [Kapfs,~]=Kfs_cal_vg(Trfs,rfsi,rnodfs,Nbe,1000,Teltfs);
    else
    
        %ask for file containing free-shape curvature in ring fixed
        %reference: first column containins angles in degrees and second 
        %one contains the freeshape curvature usual representation in 1/m
        Kapfs=dlmread('C:\Users\YangLiu\Dropbox (MIT)\Aziz Perso\papers mit\RA\Thesis\Texts\Kapfs_from_theo_pr_TC.txt'); %user has to provide this
        Tri=Kapfs(:,1)*pi/180;
        Kapfs=interp1(Tri,Kapfs(:,2),linspace(Tgg,2*pi-Tgg,1000*Nbe+1)','linear','extrap');
    end
else
    %ask for file containing force distribution giving circular shape in 
    %ring fixed reference: first column containins angles in degrees and 
    %second one contains pressure times ring height (axial width) i.e.
    %force per unit circumferential length in N/m
    PressInput=dlmread('pressure.txt'); %user has to provide this
   
    Tpo=PressInput(:,1);
    Pmo=PressInput(:,2)/Rr*Db/2;

    [~,~,Kapfs]=fscal_vg(Tpo,Pmo,Rr,Er,Izr,gap,Nbe,1000);
end

hr=hui+thrt*aui+hli+thrb*ali;
ar=arm+max(ali,aui);

 dx=1/Npe;
 x=0:dx:1;
 dxg=1/Nge;
 xg=0:dxg:1;

 N1=1-10*x.^3+15*x.^4-6*x.^5;
 N2=x-6*x.^3+8*x.^4-3*x.^5;
 N3=(x.^2)/2-3*(x.^3)/2+3*(x.^4)/2-(x.^5)/2;
 N4=10*x.^3-15*x.^4+6*x.^5;
 N5=-4*x.^3+7*x.^4-3*x.^5;
 N6=(x.^3)/2-x.^4+x.^5/2;
 N=[N1;N2;N3;N4;N5;N6];

 N1g=1-10*xg.^3+15*xg.^4-6*xg.^5;
 N2g=xg-6*xg.^3+8*xg.^4-3*xg.^5;
 N3g=(xg.^2)/2-3*(xg.^3)/2+3*(xg.^4)/2-(xg.^5)/2;
 N4g=10*xg.^3-15*xg.^4+6*xg.^5;
 N5g=-4*xg.^3+7*xg.^4-3*xg.^5;
 N6g=(xg.^3)/2-xg.^4+xg.^5/2;
 N_g=[N1g;N2g;N3g;N4g;N5g;N6g];

 Nt1=1-3*x.^2+2*x.^3;
 Nt2=x-2*x.^2+x.^3;
 Nt3=3*x.^2-2*x.^3;
 Nt4=-x.^2+x.^3;
 Nt=[Nt1;Nt2;Nt3;Nt4];

 Nt1g=1-3*xg.^2+2*xg.^3;
 Nt2g=xg-2*xg.^2+xg.^3;
 Nt3g=3*xg.^2-2*xg.^3;
 Nt4g=-xg.^2+xg.^3;
 Ntg=[Nt1g;Nt2g;Nt3g;Nt4g];
 
 dxTelt=1/min(Nge,Npe);
xTelt=0:dxTelt:1;

N1Telt=1-10*xTelt.^3+15*xTelt.^4-6*xTelt.^5;
N2Telt=Telt*(xTelt-6*xTelt.^3+8*xTelt.^4-3*xTelt.^5);
N3Telt=Telt^2*((xTelt.^2)/2-3*(xTelt.^3)/2+3*(xTelt.^4)/2-(xTelt.^5)/2);
N4Telt=10*xTelt.^3-15*xTelt.^4+6*xTelt.^5;
N5Telt=Telt*(-4*xTelt.^3+7*xTelt.^4-3*xTelt.^5);
N6Telt=Telt^2*((xTelt.^3)/2-xTelt.^4+xTelt.^5/2);

dN1Telt=-30*xTelt.^2+60*xTelt.^3-30*xTelt.^4;
dN2Telt=Telt*(1-18*xTelt.^2+32*xTelt.^3-15*xTelt.^4);
dN3Telt=Telt^2*(xTelt-9*(xTelt.^2)/2+6*xTelt.^3-5*(xTelt.^4)/2);
dN4Telt=30*xTelt.^2-60*xTelt.^3+30*xTelt.^4;
dN5Telt=Telt*(-12*xTelt.^2+28*xTelt.^3-15*xTelt.^4);
dN6Telt=Telt^2*(3*(xTelt.^2)/2-4*xTelt.^3+5*xTelt.^4/2);

d2N1Telt=-60*xTelt+180*xTelt.^2-120*xTelt.^3;
d2N2Telt=Telt*(-36*xTelt+96*xTelt.^2-60*xTelt.^3);
d2N3Telt=Telt^2*(1-9*xTelt+18*xTelt.^2-10*xTelt.^3);
d2N4Telt=60*xTelt-180*xTelt.^2+120*xTelt.^3;
d2N5Telt=Telt*(-24*xTelt+84*xTelt.^2-60*xTelt.^3);
d2N6Telt=Telt^2*(3*xTelt-12*xTelt.^2+10*xTelt.^3);

NTelt=[N1Telt;N2Telt;N3Telt;N4Telt;N5Telt;N6Telt];
dNTelt=[dN1Telt;dN2Telt;dN3Telt;dN4Telt;dN5Telt;dN6Telt];
d2NTelt=[d2N1Telt;d2N2Telt;d2N3Telt;d2N4Telt;d2N5Telt;d2N6Telt];

Nt1Telt=1-3*xTelt.^2+2*xTelt.^3;
Nt2Telt=Telt*(xTelt-2*xTelt.^2+xTelt.^3);
Nt3Telt=3*xTelt.^2-2*xTelt.^3;
Nt4Telt=Telt*(-xTelt.^2+xTelt.^3);

dNt1Telt=-6*xTelt+6*xTelt.^2;
dNt2Telt=Telt*(1-4*xTelt+3*xTelt.^2);
dNt3Telt=6*xTelt-6*xTelt.^2;
dNt4Telt=Telt*(-2*xTelt+3*xTelt.^2);

NtTelt=[Nt1Telt;Nt2Telt;Nt3Telt;Nt4Telt];
dNtTelt=[dNt1Telt;dNt2Telt;dNt3Telt;dNt4Telt];

%=====================solver======================%
prompt='Tolerance for Newton-Raphson algorithm convergence (suggested value: 1e-6): ';
Newton_tol=input(prompt); 
prompt='Maximum number of Newton-Raphson algorithm iterations (suggested value: 1000): ';
ItMax=input(prompt); 
prompt='Number of ring gap positions equally spaced between 0 and 359 degrees within bore fixed reference from thrust to thrust side to consider: ';
numTg=input(prompt);

Ug1=zeros(8*Nbnod,1);
Ug1(7:8:end-1)=alp;

F1=Finitial_fs_vg(Ug1,Er,Gr,Iyy,Izz,Jt,alp,Rr,Telt,Nbe,1000,Kapfs,dTemp,alT,hr,ar);

Kapfs=interp1(linspace(Tgg,2*pi-Tgg,1000*Nbe+1)',Kapfs,linspace(Tgg,2*pi-Tgg,Nbe*min(Nge,Npe)+1)','linear','extrap');

kapfsz=Kapfs*sin(alp);

inputsl4.ring.F1=F1;

[Kg_nd]=Kg_ring_curved_nd_vg(Er,Gr,Iyy,Jt,Rr,Telt,Nbe);

inputsl4.FEM.Kg_nd=Kg_nd;

UgIC_nd=zeros(8*Nbnod,1);

 Nplot=min(Ng,Nr);
 res_Tg=linspace(0,359,numTg)';
 
 res_ybr=zeros(Nr+1,numTg);
 res_hmin=zeros(Nr+1,numTg);
 res_zr=zeros(Nr+1,numTg);
 res_yr=zeros(Nr+1,numTg);
 res_alr=zeros(Nr+1,numTg);
 res_z0=zeros(Nr+1,numTg);
 res_uoc=zeros(Nr+1,numTg);
 res_uic=zeros(Nr+1,numTg);
 res_loc=zeros(Nr+1,numTg);
 res_lic=zeros(Nr+1,numTg);
 
 res_fy=zeros(Nplot+1,numTg);
 res_fl=zeros(Nplot+1,numTg);
 res_fz=zeros(Nplot+1,numTg);
 res_m=zeros(Nplot+1,numTg);
 res_Mfinal=zeros(Nplot+1,numTg);
 res_S_uo=zeros(Nplot+1,numTg);
 res_S_lo=zeros(Nplot+1,numTg);
 res_S_ui=zeros(Nplot+1,numTg);
 res_S_li=zeros(Nplot+1,numTg);

 util_data=[Nplot Nr Db gap IsBDist];

 for jframe=1:numTg

Tg=res_Tg(jframe)/180*pi;   
Tnod=(Tgg:(Telt):2*pi-Tgg)+Tg;
Vnod=[Rr*cos(Tnod),Rr*sin(Tnod)];
shiftstep=round(Tg/2/pi*Ng);
idg2r=circshift(id,-shiftstep);
idr2g=circshift(id,shiftstep);

shiftstepr=round(Tg/2/pi*Nr);
shiftstepg=round(Tg/2/pi*Ng);
mu_oil=circshift(mu_oilglob,[0 -shiftstepg]);
mu_oil=interp1(anglesg,mu_oil,angle_thetag,'linear','extrap');

Pi=circshift(Piglob,[0 -shiftstepg]);
Pi=interp1(anglesg,Pi,angle_thetag,'linear','extrap');
Pd=circshift(Pdglob,[0 -shiftstepg]);
Pd=interp1(anglesg,Pd,angle_thetag,'linear','extrap');
Pu=circshift(Puglob,[0 -shiftstepg]);
Pu=interp1(anglesg,Pu,angle_thetag,'linear','extrap');

hrls=circshift(hrlsglob,[0 -shiftstepr]);
hrls=interp1(anglesr,hrls,angle_theta,'linear','extrap');
hot=circshift(hotglob,[0 -shiftstepg]);
hot=interp1(anglesg,hot,angle_thetag,'linear','extrap');
hob=circshift(hobglob,[0 -shiftstepg]);
hob=interp1(anglesg,hob,angle_thetag,'linear','extrap');

templ=circshift(templglob,[0 -shiftstepr]);
templ=interp1(anglesr,templ,angle_theta,'linear','extrap');
isfuel=circshift(isfuelglob,[0 -shiftstepr]);
isfuel=interp1(anglesr,isfuel,angle_theta,'linear','extrap');
isfuel=logical(isfuel);

inputsl4.ring.idr2g=idr2g;
inputsl4.ring.idg2r=idg2r;
inputsl4.ring.Tg=Tg;
inputsl4.FEM.Tnod=Tnod;
inputsl4.FEM.Vnod=Vnod;

inputsl4.oil.mu_oil=mu_oil;

inputsl4.engine.Pu=Pu;
inputsl4.engine.Pd=Pd;
inputsl4.engine.Pi=Pi;

inputsl4.lubrication.hrls=hrls;
inputsl4.lubrication.hot=hot;
inputsl4.lubrication.hob=hob;
inputsl4.lubrication.isfuel=isfuel;

inputsl4.bore.templ=templ;

[ybn]=bore_dist_static_vg(inputsl4);
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

UgIC_nd(1:8:8*Nbnod-7)=ybn(id1)+(Db/2-Rr-arm-sigmap)/uref;
UgIC_nd(2:8:8*Nbnod-6)=ybn(id2);
UgIC_nd(3:8:8*Nbnod-5)=ybn(id3);

Pua=mean(Pu);
Pia=mean(Pi);
Pda=mean(Pd);
dz_g=(Rr*cos(Tnod')-off)*beta_p;
[zgnu,zgnl]=groove_dist_static_vg(inputsl4);

zgnu1=zgnu(id1);
zgnu2=zgnu(id2);
zgnu3=zgnu(id3);

zgnu(id1)=interp1(Tnodd,zgnu1,Tnod');
zgnu(id2)=interp1(Tnodd,zgnu2,Tnod');
zgnu(id3)=interp1(Tnodd,zgnu3,Tnod');

zgnl1=zgnl(id1);
zgnl2=zgnl(id2);
zgnl3=zgnl(id3);

zgnl(id1)=interp1(Tnodd,zgnl1,Tnod');
zgnl(id2)=interp1(Tnodd,zgnl2,Tnod');
zgnl(id3)=interp1(Tnodd,zgnl3,Tnod');

hg_ig=(omega-((-F1(4:8:end-4)*2/Le/lfu+Fi/2/pi/Rr/lfu)/Pk_rg+Pda-Pua).^(1/z))*sigmagt.*((F1(4:8:end-4)*2/Le/lfu-Fi/2/pi/Rr/lfu+Pua-Pda)<0)...
     +(omega-((F1(4:8:end-4)*2/Le/lfl-Fi/2/pi/Rr/lfl+Pua-Pda)/Pk_rg).^(1/z))*sigmagb.*((F1(4:8:end-4)*2/Le/lfl-Fi/2/pi/Rr/lfl+Pua-Pda)>=0);
UgIC_nd(4:8:8*Nbnod-4)=((gcl/2-hg_ig+dz_g)/href+zgnu(id1)).*((F1(4:8:end-4)*2/Le/lfu-Fi/2/pi/Rr/lfu+Pua-Pda)<0)...
     +((hg_ig-gcl/2+dz_g)/href+zgnl(id1)).*((F1(4:8:end-4)*2/Le/lfl-Fi/2/pi/Rr/lfl+Pua-Pda)>=0);
UgIC_nd(5:8:8*Nbnod-3)=(-Rr*sin(Tnod')*beta_p/href*Telt+zgnu(id2)).*((F1(4:8:end-4)*2/Le/lfu-Fi/2/pi/Rr/lfu+Pua-Pda)<0)...
     +(-Rr*sin(Tnod')*beta_p/href*Telt+zgnl(id2)).*((F1(4:8:end-4)*2/Le/lfl-Fi/2/pi/Rr/lfl+Pua-Pda)>=0);
UgIC_nd(6:8:8*Nbnod-2)=(-Rr*cos(Tnod')*beta_p/href*Telt^2+zgnu(id3)).*((F1(4:8:end-4)*2/Le/lfu-Fi/2/pi/Rr/lfu+Pua-Pda)<0)...
     +(-Rr*cos(Tnod')*beta_p/href*Telt^2+zgnl(id3)).*((F1(4:8:end-4)*2/Le/lfl-Fi/2/pi/Rr/lfl+Pua-Pda)>=0);

if Pua>(Pia+Pda)
    beta=beta_p*cos(Tnod')+tilt_thu+1e-6;
    UgIC_nd(7:8:8*Nbnod-1)=beta/href*Rr-sign(F1(7))*1e-6*sin((Tgg:Telt:2*pi-Tgg)/2)'/href*Rr;
    UgIC_nd(8:8:8*Nbnod)=-beta_p*sin(Tnod')/href*Rr*Telt-sign(F1(7))*1e-6*cos((Tgg:Telt:2*pi-Tgg)'/2)/href*Rr*Telt;

else
    UgIC_nd(7:8:8*Nbnod-1)=-sign(F1(7))*1e-4*sin((Tgg:Telt:2*pi-Tgg)/2)/href*Rr;
    UgIC_nd(8:8:8*Nbnod)=-sign(F1(7))*1e-4*cos((Tgg:Telt:2*pi-Tgg)/2)/href*Rr*Telt;
end

Ug_nd=UgIC_nd;

judge1=1;
judge2=1;
[NLS_0,yr,zr,alr,hmin,~,z0,~,~,~]=Static_Error_vg(Ug_nd,F1,inputs);
NLS_k=NLS_0;
j=0;
while j<=ItMax&&(judge1||judge2)
    j=j+1;
    [Ja_stac]=Jacobian_Estac_vg(Ug_nd,inputsl4);
    Nstep=-Ja_stac\NLS_k;
    Ugo_nd=Ug_nd;
    Erro=norm(NLS_k);
    itGB=0; dErr=-1;
    
    if j<4
        while (dErr<0)&&(itGB<20)
            if itGB>0
                Nstep=Nstep/4;
            end
            Ug_nd=Ugo_nd+Nstep; 
            [NLS_k,yr,zr,alr,hmin,~,z0,~,~,~]=Static_Error_vg(Ug_nd,F1,inputs);
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
            [NLS_k,yr,zr,alr,hmin,~,z0,~,~,~]=Static_Error_vg(Ug_nd,F1,inputs);
            Err=norm(NLS_k); 
            dErr=Erro-Err;
            itGB=itGB+1;
        end
    end

    judge1=norm(Nstep)>Newton_tol;
    judge2=norm(NLS_k)>Newton_tol;
    
end

hmin=max(0,hmin);

 [ybn]=bore_dist_static_vg(inputsl4);
 Tnodd=(0:2*pi/Nbe:2*pi)+Tg;
 
id1=1:3:3*Nbnod-2;
id2=2:3:3*Nbnod-1;
id3=3:3:3*Nbnod;

ybn1=ybn(id1);
ybn2=ybn(id2);
ybn3=ybn(id3);

ybn(id1)=interp1(Tnodd,ybn1,Tnod');
ybn(id2)=interp1(Tnodd,ybn2,Tnod');
ybn(id3)=interp1(Tnodd,ybn3,Tnod');

 [zgnu,zgnl]=groove_dist_static_vg(inputsl4);

 zgnu1=zgnu(id1);
zgnu2=zgnu(id2);
zgnu3=zgnu(id3);

zgnu(id1)=interp1(Tnodd,zgnu1,Tnod');
zgnu(id2)=interp1(Tnodd,zgnu2,Tnod');
zgnu(id3)=interp1(Tnodd,zgnu3,Tnod');

zgnl1=zgnl(id1);
zgnl2=zgnl(id2);
zgnl3=zgnl(id3);

zgnl(id1)=interp1(Tnodd,zgnl1,Tnod');
zgnl(id2)=interp1(Tnodd,zgnl2,Tnod');
zgnl(id3)=interp1(Tnodd,zgnl3,Tnod');

 ybr=zeros(1,Nr+1);
 yrr=zeros(1,Nr+1);
 uoc=zeros(1,Nr+1);
 uic=zeros(1,Nr+1);
 loc=zeros(1,Nr+1);
 lic=zeros(1,Nr+1);
  
 for j=1:Nbe
     ind=1+(j-1)*Npe:1+j*Npe;
     
     Uey=Ug_nd([8*(j-1)+1:8*(j-1)+3,8*j+1:8*j+3]);
     Uez=Ug_nd([8*(j-1)+4:8*(j-1)+6,8*j+4:8*j+6]);
     Uet=Ug_nd([8*j-1:8*j,8*j+7:8*j+8]);
     
     yrl=(Uey')*N*uref;
     zrl=(Uez')*N*href;
     alrl=(Uet')*Nt*href/Rr;
     beta=beta_p*cos((Tnod(j)+Telt*x));
     beta_thu=tilt_thu*ones(size(x)); 
     beta_thl=tilt_thl*ones(size(x));
     dz_g=(Rr*cos((Tnod(j)+Telt*x))-off)*beta_p;
     
     id=3*(j-1)+1:3*(j+1);
     yb=(ybn(id))'*N*uref;
     
     ybr(ind)=yb;
     yrr(ind)=yrl;
     
     zgu=(zgnu(id))'*N*href;
     zgl=(zgnl(id))'*N*href; 
     
     Nby=50;

     yu=linspace(y1u,y2u,Nby+1);
     yl=linspace(y1l,y2l,Nby+1);
     
     z_2u=repmat(dz_g+zgu-zrl,Nby+1,1);
     al_2u=repmat(beta+beta_thu-alrl,Nby+1,1);
     z_2l=repmat(dz_g+zgl-zrl,Nby+1,1);
     al_2l=repmat(beta+beta_thl-alrl,Nby+1,1);
 
     
     yyu=repmat(yu',1,Npe+1);
     yyl=repmat(yl',1,Npe+1);
     
     hgu=gcl/2+z_2u+al_2u.*yyu;
     hgl=gcl/2-z_2l-al_2l.*yyl;
     uoc(ind)=hgu(end,:);
     uic(ind)=hgu(1,:);
     loc(ind)=hgl(end,:);
     lic(ind)=hgl(1,:);
 
 end
 
 theer=linspace(Tgg*180/pi,360-Tgg*180/pi,Nr+1)';

 uoc=max(0,uoc);
 uic=max(0,uic);
 loc=max(0,loc);
 lic=max(0,lic);
 
 res_theer=theer;
 res_hmin(:,jframe)=hmin'*1e6;
 res_zr(:,jframe)=zr'*1e6;
 res_yr(:,jframe)=yr'*1e6;
 res_alr(:,jframe)=alr'/pi*180;
 res_z0(:,jframe)=z0'*1e6;
 res_uoc(:,jframe)=uoc'*1e6;
 res_uic(:,jframe)=uic'*1e6;
 res_loc(:,jframe)=loc'*1e6;
 res_lic(:,jframe)=lic'*1e6;
 temp=inputsl4.bore.templ;
 dt=1;
 sigmagt=inputsl4.piston.sigmag1t;
 sigmagb=inputsl4.piston.sigmag1b;
 Ph_OCR=inputsl4.liner.Ph;
 fc=inputsl4.contact.fc_dry;
 isfuela=inputsl4.lubrication.isfuel;
 temp_liner=temp;
 
 isfuel=isfuela;
 hrls=hrls.*(hrls>=(minoil))+minoil*(hrls<(minoil));
 hrl=2*hrls;
 Kp=cp1+cp2*hrl/sigmap;
 Kf=cf1+cf2*hrl/sigmap;

fcyglob=zeros(1,Nr+1);
fczglob=zeros(1,Nr+1);
fgyglob=zeros(1,Ng+1);
fgzglob=zeros(1,Ng+1);
mcglob=zeros(1,Nr+1);
mgglob=zeros(1,Ng+1);
fhlglob=zeros(1,Nr+1);
fclglob=zeros(1,Nr+1);
for j=1:Nbe
    
    ind=(Npe*(j-1)+1):(Npe*j+1);
    indg=(Nge*(j-1)+1):(Nge*j)+1;
    Uey=Ug_nd([8*(j-1)+1:8*(j-1)+3,8*j+1:8*j+3]);
    Uez=Ug_nd([8*(j-1)+4:8*(j-1)+6,8*j+4:8*j+6]);
    Uet=Ug_nd([8*j-1:8*j,8*j+7:8*j+8]);
    
    yre=(Uey')*N*uref;
    zre=(Uez')*N*href;
    zrg=(Uez')*N_g*href;
    zrgm1=zrg;
    alre=(Uet')*Nt*href/Rr;
    alrg=(Uet')*Ntg*href/Rr;
    alrgm1=alrg;
    beta=beta_p*cos((Tnod(j)+Telt*x));
    betag=beta_p*cos((Tnod(j)+Telt*xg));
    betagm1=betag;
    
    beta_thu=tilt_thu*ones(size(x));
    beta_thug=tilt_thu*ones(size(xg));
    
    beta_thl=tilt_thl*ones(size(x));
    beta_thlg=tilt_thl*ones(size(xg));
    dz_g=(Rr*cos((Tnod(j)+Telt*x))-off)*beta_p;
    dz_gg=(Rr*cos((Tnod(j)+Telt*xg))-off)*beta_p;
    dz_ggm1=dz_gg;
    
    id=3*(j-1)+1:3*(j+1);
    yb=(ybn(id))'*N*uref;
    
    zgu=(zgnu(id))'*N*href;  
    zgl=(zgnl(id))'*N*href; 
    zgug=(zgnu(id))'*N_g*href;  
    zglg=(zgnl(id))'*N_g*href;  
    
    templ=temp_liner(ind);
    hrle=hrl(ind);
    Kpe=Kp(ind);
    Kfe=Kf(ind);
    
    hote=hot(indg);
    hobe=hob(indg);
    mu_oile=mu_oil(indg);
    Pue=Pu(indg);
    Pde=Pd(indg);
    Pie=Pi(indg);
    isfuele=isfuel(ind);
    
    [fcgu,mcgu]=ring_groove_contact_vg(gcl,dz_g,zgu,zre,alre,beta,(beta_thu+thgt-thrt),y1u,y2u,sigmagt,z,Pk_rg,omega,1);
    [fcgl,mcgl]=ring_groove_contact_vg(gcl,dz_g,zgl,zre,alre,beta,(beta_thl-thgb+thrb),y1l,y2l,sigmagb,z,Pk_rg,omega,-1);
    [fcl,mcl,~]=ring_liner_contact_vg(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yre,alre,PR,sigmap,z,Pk_rl,omega);
    [~,~,fhl,ffl,mhl]=ring_liner_hydro_vg(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yre,alre,hrle,KOCR,Kpe,Kfe,Ph_OCR,ap,F0,b,c,d,e,rw,muUot,templ,isfuele,Vp,inputsl4.oil,sigmap,oilthreshold,viscosityfactor);
    [frgas,fzgas,mgas]=ring_gas_exg_vg(Pie,Pue,Pde,arm,aui,lfu,ali,lfl,hui,hli,huo,hlo,rbn,rb1,rb2,a11,a12,a21,a22,alrg);  
    [fgu,mgu]=ring_tgroove_og_vg(Pie,Pue,gcl,dz_gg,dz_ggm1,zgug,zrg,zrgm1,alrg,alrgm1,betag,betagm1,(beta_thug+thgt-thrt),y1u,y2u,hote,mu_oile,dt);
    [fgl,mgl]=ring_bgroove_og_vg(Pie,Pde,gcl,dz_gg,dz_ggm1,zglg,zrg,zrgm1,alrg,alrgm1,betag,betagm1,(beta_thlg-thgb+thrb),y1l,y2l,hobe,mu_oile,dt);   
    
    directionf=-sign(Vp);
    ffc=(-fcl)*fc*directionf;
    
    fcz=ffc+fcgu+fcgl+ffl;
    fgz=fgu+fgl+fzgas;
    
    fczglob(ind)=fcz;
    fgzglob(indg)=fgz;

    fcy=-fcgu*thrt+fcgl*thrb+fcl+fhl;
    fgy=-fgu*thrt+fgl*thrb+frgas;
    fhlglob(ind)=fhl;
    fclglob(ind)=fcl;
    fcyglob(ind)=fcy;
    fgyglob(indg)=fgy;

    mca=mcgu+fcgu*hui*thrt+mcgl+fcgl*hli*thrb+mcl+mhl+ffc*arm;
    mga=mgu+fgu*hui*thrt+mgl+fgl*hli*thrb+mgas;
    
    mcglob(ind)=mca;
    mgglob(indg)=mga;

end

 theeg=linspace(Tgg*180/pi,360-Tgg*180/pi,Ng+1);
 theer=linspace(Tgg*180/pi,360-Tgg*180/pi,Nr+1);
 thee=linspace(Tgg*180/pi,360-Tgg*180/pi,Nplot+1)';
 fgyglobplot=interp1(theeg,fgyglob,thee);
 fcyglob=interp1(theer,fcyglob,thee);
 fczglob=interp1(theer,fczglob,thee);
 mcglob=interp1(theer,mcglob,thee);
 fhlglob=interp1(theer,fhlglob,thee);
 fclglob=interp1(theer,fclglob,thee);
 fgzglobplot=interp1(theeg,fgzglob,thee);
 mgglobplot=interp1(theeg,mgglob,thee);
 fy=fcyglob+fgyglobplot;
 fz=fczglob+fgzglobplot;
 m=mcglob+mgglobplot;
 fl=fhlglob+fclglob;

 res_thee=thee;
res_fy(:,jframe)=-fy;
res_fl(:,jframe)=-fl;
res_fz(:,jframe)=fz;
res_m(:,jframe)=m;

Uov=zeros(1,3*Nbnod);
Uezn=zeros(1,3*Nbnod);
Uean=zeros(1,2*Nbnod);

for j=1:Nbe+1
     Uov(3*(j-1)+1:3*(j-1)+3)=uref*Ug_nd(8*(j-1)+1:8*(j-1)+3);
     Uezn(3*(j-1)+1:3*(j-1)+3)=href*Ug_nd(8*(j-1)+4:8*(j-1)+6);
     Uean(3*(j-1)+1:3*(j-1)+2)=href/Rr*Ug_nd(8*(j-1)+7:8*(j-1)+8);
end
Uov=Uov';
Uezn=Uezn';
Uean=Uean';

ov=zeros(1,Nplot+1);
dov=zeros(1,Nplot+1);
d2ov=zeros(1,Nplot+1);

dz=zeros(1,Nplot+1);
d2z=zeros(1,Nplot+1);

alpha=zeros(1,Nplot+1);
dalpha=zeros(1,Nplot+1);

for i=1:Nbe-1
    
    Uey=Uov(3*i-2:3*i+3,1);
    Uez=Uezn(3*i-2:3*i+3);
    Uea=Uean(3*i-2:3*i+1);
    
    ov(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uey')*NTelt(:,1:end-1);
    dov(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uey')*dNTelt(:,1:end-1)/Telt;
    d2ov(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uey')*d2NTelt(:,1:end-1)/Telt^2;
    
    dz(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uez')*dNTelt(:,1:end-1)/Telt;
    d2z(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uez')*d2NTelt(:,1:end-1)/Telt^2;
    
    alpha(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uea')*NtTelt(:,1:end-1);
    dalpha(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uea')*dNtTelt(:,1:end-1)/Telt;
    
end

Uey=Uov(end-5:end);
ov(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uey')*NTelt;
dov(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uey')*dNTelt/Telt;
d2ov(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uey')*d2NTelt/Telt^2;

Uez=Uezn(end-5:end);
dz(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uez')*dNTelt/Telt;
d2z(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uez')*d2NTelt/Telt^2;

Uea=Uean(end-3:end);
alpha(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uea')*NtTelt;
dalpha(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uea')*dNtTelt/Telt;

kzz=(1/Rr-(ov+d2ov)/Rr^2).*sin(alpha+d2z/Rr);
kzz=kzz';

ov=ov';
ov2=ov(Nplot/2+1:end);
ov2=flipud(ov2);
ov=ov(1:Nplot/2+1);
ov_uo=ov+auo;
ov_uo2=ov2+auo;
ov_uo2=flipud(ov_uo2);
ov_lo=ov+alo;
ov_lo2=ov2+alo;
ov_lo2=flipud(ov_lo2);
ov_ui=ov-aui;
ov_ui2=ov2-aui;
ov_ui2=flipud(ov_ui2);
ov_li=ov-ali;
ov_li2=ov2-ali;
ov_li2=flipud(ov_li2);
dov=dov';
dov2=dov(Nplot/2+1:end);
dov2=flipud(dov2);
dov=dov(1:Nplot/2+1);
Tp=linspace(Tgg,pi,Nplot/2+1)';
Tgloc=Tgg;

Pm=-fy;
Pm2=Pm(Nplot/2+1:end);
Pm=Pm(1:Nplot/2+1);

Mfinal=zeros(length(Tp),1);
Mfinal2=zeros(length(Tp),1);
S_uo=zeros(length(Tp),1);
S_uo2=zeros(length(Tp),1);
S_lo=zeros(length(Tp),1);
S_lo2=zeros(length(Tp),1);
S_ui=zeros(length(Tp),1);
S_ui2=zeros(length(Tp),1);
S_li=zeros(length(Tp),1);
S_li2=zeros(length(Tp),1);

for i=2:length(Tp)
    n=1000;
    Tha=linspace(Tgloc,Tp(i),n)';

    uov=interp1(Tp,ov,Tha);
    uov_uo=interp1(Tp,ov_uo,Tha);
    uov_lo=interp1(Tp,ov_lo,Tha);
    uov_ui=interp1(Tp,ov_ui,Tha);
    uov_li=interp1(Tp,ov_li,Tha);

    uov2=interp1(Tp,ov2,Tha);
    uov_uo2=interp1(Tp,ov_uo2,Tha);
    uov_lo2=interp1(Tp,ov_lo2,Tha);
    uov_ui2=interp1(Tp,ov_ui2,Tha);
    uov_li2=interp1(Tp,ov_li2,Tha);

    udov=interp1(Tp,dov,Tha);
    udov2=interp1(Tp,dov2,Tha);
    da=(Tp(i)-Tgloc)/n;
    upov=interp1(Tp,Pm,Tha);
    upov2=interp1(Tp,Pm2,Tha);

    df=abs((udov.*cos(Tha)-(Rr+uov).*sin(Tha)).*(Rr+ov(i)).*cos(Tp(i))+(udov.*sin(Tha)+(Rr+uov).*cos(Tha)).*(Rr+ov(i)).*sin(Tp(i))-(Rr+uov).*udov)./sqrt((udov.*cos(Tha)-(Rr+uov).*sin(Tha)).^2+(udov.*sin(Tha)+(Rr+uov).*cos(Tha)).^2);
    df_uo=abs((udov.*cos(Tha)-(Rr+uov_uo).*sin(Tha)).*(Rr+ov(i)).*cos(Tp(i))+(udov.*sin(Tha)+(Rr+uov_uo).*cos(Tha)).*(Rr+ov(i)).*sin(Tp(i))-(Rr+uov_uo).*udov)./sqrt((udov.*cos(Tha)-(Rr+uov_uo).*sin(Tha)).^2+(udov.*sin(Tha)+(Rr+uov_uo).*cos(Tha)).^2);
    df_lo=abs((udov.*cos(Tha)-(Rr+uov_lo).*sin(Tha)).*(Rr+ov(i)).*cos(Tp(i))+(udov.*sin(Tha)+(Rr+uov_lo).*cos(Tha)).*(Rr+ov(i)).*sin(Tp(i))-(Rr+uov_lo).*udov)./sqrt((udov.*cos(Tha)-(Rr+uov_lo).*sin(Tha)).^2+(udov.*sin(Tha)+(Rr+uov_lo).*cos(Tha)).^2);
    df_ui=abs((udov.*cos(Tha)-(Rr+uov_ui).*sin(Tha)).*(Rr+ov(i)).*cos(Tp(i))+(udov.*sin(Tha)+(Rr+uov_ui).*cos(Tha)).*(Rr+ov(i)).*sin(Tp(i))-(Rr+uov_ui).*udov)./sqrt((udov.*cos(Tha)-(Rr+uov_ui).*sin(Tha)).^2+(udov.*sin(Tha)+(Rr+uov_ui).*cos(Tha)).^2);
    df_li=abs((udov.*cos(Tha)-(Rr+uov_li).*sin(Tha)).*(Rr+ov(i)).*cos(Tp(i))+(udov.*sin(Tha)+(Rr+uov_li).*cos(Tha)).*(Rr+ov(i)).*sin(Tp(i))-(Rr+uov_li).*udov)./sqrt((udov.*cos(Tha)-(Rr+uov_li).*sin(Tha)).^2+(udov.*sin(Tha)+(Rr+uov_li).*cos(Tha)).^2);
    
    df2=abs((udov2.*cos(Tha)-(Rr+uov2).*sin(Tha)).*(Rr+ov2(i)).*cos(Tp(i))+(udov2.*sin(Tha)+(Rr+uov2).*cos(Tha)).*(Rr+ov2(i)).*sin(Tp(i))-(Rr+uov2).*udov2)./sqrt((udov2.*cos(Tha)-(Rr+uov2).*sin(Tha)).^2+(udov2.*sin(Tha)+(Rr+uov2).*cos(Tha)).^2);
    df_uo2=abs((udov2.*cos(Tha)-(Rr+uov_uo2).*sin(Tha)).*(Rr+ov2(i)).*cos(Tp(i))+(udov2.*sin(Tha)+(Rr+uov_uo2).*cos(Tha)).*(Rr+ov2(i)).*sin(Tp(i))-(Rr+uov_uo2).*udov2)./sqrt((udov2.*cos(Tha)-(Rr+uov_uo2).*sin(Tha)).^2+(udov2.*sin(Tha)+(Rr+uov_uo2).*cos(Tha)).^2);
    df_lo2=abs((udov2.*cos(Tha)-(Rr+uov_lo2).*sin(Tha)).*(Rr+ov2(i)).*cos(Tp(i))+(udov2.*sin(Tha)+(Rr+uov_lo2).*cos(Tha)).*(Rr+ov2(i)).*sin(Tp(i))-(Rr+uov_lo2).*udov2)./sqrt((udov2.*cos(Tha)-(Rr+uov_lo2).*sin(Tha)).^2+(udov2.*sin(Tha)+(Rr+uov_lo2).*cos(Tha)).^2);
    df_ui2=abs((udov2.*cos(Tha)-(Rr+uov_ui2).*sin(Tha)).*(Rr+ov2(i)).*cos(Tp(i))+(udov2.*sin(Tha)+(Rr+uov_ui2).*cos(Tha)).*(Rr+ov2(i)).*sin(Tp(i))-(Rr+uov_ui2).*udov2)./sqrt((udov2.*cos(Tha)-(Rr+uov_ui2).*sin(Tha)).^2+(udov2.*sin(Tha)+(Rr+uov_ui2).*cos(Tha)).^2);
    df_li2=abs((udov2.*cos(Tha)-(Rr+uov_li2).*sin(Tha)).*(Rr+ov2(i)).*cos(Tp(i))+(udov2.*sin(Tha)+(Rr+uov_li2).*cos(Tha)).*(Rr+ov2(i)).*sin(Tp(i))-(Rr+uov_li2).*udov2)./sqrt((udov2.*cos(Tha)-(Rr+uov_li2).*sin(Tha)).^2+(udov2.*sin(Tha)+(Rr+uov_li2).*cos(Tha)).^2);
    
    dMe=upov.*(Rr+uov).*df*da;
    dMe_uo=upov.*(Rr+uov_uo).*df_uo*da;
    dMe_lo=upov.*(Rr+uov_lo).*df_lo*da;
    dMe_ui=upov.*(Rr+uov_ui).*df_ui*da;
    dMe_li=upov.*(Rr+uov_li).*df_li*da;
    
    dMe2=upov2.*(Rr+uov2).*df2*da;
    dMe_uo2=upov2.*(Rr+uov_uo2).*df_uo2*da;
    dMe_lo2=upov2.*(Rr+uov_lo2).*df_lo2*da;
    dMe_ui2=upov2.*(Rr+uov_ui2).*df_ui2*da;
    dMe_li2=upov2.*(Rr+uov_li2).*df_li2*da;

    Mfinal(i)=sum(dMe(1:end-1)+dMe(2:end))/2;
    S_uo(i)=auo/Izr*sum(dMe_uo(1:end-1)+dMe_uo(2:end))/2;
    S_lo(i)=alo/Izr*sum(dMe_lo(1:end-1)+dMe_lo(2:end))/2;
    S_ui(i)=-aui/Izr*sum(dMe_ui(1:end-1)+dMe_ui(2:end))/2;
    S_li(i)=-ali/Izr*sum(dMe_li(1:end-1)+dMe_li(2:end))/2;
    
    Mfinal2(i)=sum(dMe2(1:end-1)+dMe2(2:end))/2;
    S_uo2(i)=auo/Izr*sum(dMe_uo2(1:end-1)+dMe_uo2(2:end))/2;
    S_lo2(i)=alo/Izr*sum(dMe_lo2(1:end-1)+dMe_lo2(2:end))/2;
    S_ui2(i)=-aui/Izr*sum(dMe_ui2(1:end-1)+dMe_ui2(2:end))/2;
    S_li2(i)=-ali/Izr*sum(dMe_li2(1:end-1)+dMe_li2(2:end))/2;

end

 Mfinal=[Mfinal(1:end-1);flipud(Mfinal2)];
 S_uo=[S_uo(1:end-1);flipud(S_uo2)];
 ruo=sqrt(auo^2+huo^2);
 S_uo=S_uo-abs(Er*Iyy*(kzz-kapfsz))*huo/Iyr.*(fz>=0)+abs(Er*Iyy*(kzz-kapfsz))*huo/Iyr.*(fz>=0)+ruo*Gr*(dz'/Rr^2-dalpha'/Rr);
 S_uo=S_uo*1e-6;
 
 S_lo=[S_lo(1:end-1);flipud(S_lo2)];
 rlo=sqrt(alo^2+hlo^2);
 S_lo=S_lo+abs(Er*Iyy*(kzz-kapfsz))*hlo/Iyr.*(fz>=0)-abs(Er*Iyy*(kzz-kapfsz))*hlo/Iyr.*(fz>=0)+rlo*Gr*(dz'/Rr^2-dalpha'/Rr);
 S_lo=S_lo*1e-6;
 
 S_ui=[S_ui(1:end-1);flipud(S_ui2)];
 rui=sqrt(aui^2+hui^2);
 S_ui=S_ui-abs(Er*Iyy*(kzz-kapfsz))*hui/Iyr.*(fz>=0)+abs(Er*Iyy*(kzz-kapfsz))*hui/Iyr.*(fz>=0)+rui*Gr*(dz'/Rr^2-dalpha'/Rr);
 S_ui=S_ui*1e-6;
 
 S_li=[S_li(1:end-1);flipud(S_li2)];
 rli=sqrt(ali^2+hli^2);
 S_li=S_li+abs(Er*Iyy*(kzz-kapfsz))*hli/Iyr.*(fz>=0)-abs(Er*Iyy*(kzz-kapfsz))*hli/Iyr.*(fz>=0)+rli*Gr*(dz'/Rr^2-dalpha'/Rr);
 S_li=S_li*1e-6;
 
 res_Mfinal(:,jframe)=Mfinal;
 res_S_uo(:,jframe)=S_uo;
 res_S_lo(:,jframe)=S_lo;
 res_S_ui(:,jframe)=S_ui;
 res_S_li(:,jframe)=S_li;
 
 [ybn]=bore_dist_static_vg(inputsl4);
 for j=1:Nbe
     ind=1+(j-1)*Npe:1+j*Npe;
     id=3*(j-1)+1:3*(j+1);
     yb=(ybn(id))'*N*uref;
     
     ybr(ind)=yb;
 end
 res_ybr(:,jframe)=ybr';
 
 end
 
end