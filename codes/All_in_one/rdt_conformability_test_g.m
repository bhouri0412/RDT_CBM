function rdt_conformability_test_g()

%this code performs the computation related to the conformability study for on specific ring gap position
%the results are given in the variables mat_res_1 which is a matrix whose
%columns contain from left to right: angles in circumference direction
%(degrees), rin-liner minimum clearance (microns), axial lift (microns),
%radial coordinate with respect to ring nominal radius (microns), twist
%angles (degrees), minimum clearance point axial location % centroid
%(microns), upper OD clearance (microns), upper ID clearance (microns),
%lower OD clearance (microns), lower ID clearance (microns)

% the rest of the results are given in the variable mat_res_2 which is a matrix whose
%columns contain from left to right: angles in circumference direction
%(degrees), radial force (N/m), liner force (N/m), axial force (N/m),
%twist moment (N), bending moment (Nm), stress at upper OD (Nm), stress at
%upper ID (Nm), stress at lower OD (Nm), stress at lower ID (Nm)

close all

addpath(genpath(pwd))

%=============features================================%
IsBDist=0; %1 to take into account bore thermal distortion, 0 otherwise
IsPTilt=0; %1 to take into account piston dynamic tilt, 0 otherwise
IsGDef=0; %1 to take into account groove thermal distortion, 0 otherwise
IsFreeshape=0;  %1 to provide free shape or its curvature and 0 to provide radial pressure distribution for the circular shape

mag1=3000; %Magnification coefficient for radial plot of bore distortion and ring liner clearance

numforce=200; %Number of radial forces to represent (suggested value: 200)
fstop=20000; %Maximum acceptable value for forces to be represented (N/m) (suggested value: 20000)

mag2=mag1;

inputs.features.IsBDist=IsBDist;
inputs.features.IsPTilt=IsPTilt;
inputs.features.IsGDef=IsGDef;
inputs.features.IsFreeshape=IsFreeshape;

%============================ring input=======================%
Db=95.25e-3;  %bore diameter, m
auo=1.939644012944984e-3;   %upper OD width, m
aui=2.060355987055016e-3;   %upper ID width, m
alo=1.939644012944984e-3;   %lower OD width, m
ali=2.060355987055016e-3;   %lower ID width, m
huo=1e-3;   %upper OD height, m
hui=1e-3;  %upper ID height, m
hlo=1e-3;  %lower OD height, m
hli=1e-3;  %lower ID height, m
thrt=0;           %ring upper flank angle, m
thrb=0;           %ring lower flank angle, m
 
rb1=0.7e-3;    %the lower edge width, m
rb2=0.7e-3;     %the upper edge width, m
rbn=-0.00227088e-3;  %running face minimum point axial location, m
a10=0; 
a11=0; %Linear cofficient for lower edge shape factor:
a12=50;        %the lower edge shape factor, 1/m
a20=0;
a21=0; %Linear cofficient for upper edge shape factor:
a22=50;       %the upper edge shape factor, 1/m

arm=2.139644012944984e-3; %running face minimum point width, m
gap=0.334e-3; %ring gap size when closed inside bore, m

Ac=8.22e-6;   %area of cross-section, m^2
Izz=11.667782955771303e-12; %principal moment of inertial in plane, m^4
Iyy=2.714600000000000e-12; %principal moment of inertial out of plane, m^4
Iyr=2.714600000000000e-12; %moment of inertial out of plane, m^4
Izr=11.667782955771303e-12; %moment of inertial in plane, m^4
Ip=Izz+Iyy;   
Jt=7.573654782911500e-12;   %torsional factor, m^4
alp=0/180*pi;  %principal angle, rad

rhor=7820;       %ring density, kg/m^3
Er=100e9;        %ring Young's modulus, Pa
nur=0.3;         %ring Poisson Ratio
Gr=Er/(2*(1+nur));
alT=1.1e-5;       %thermal expansion coefficient, 1/K

Tg=0/180*pi;    %ring gap position, rad

rTemp=25;        %ring temperature, C

Rro=Db/2-arm;
Rr=Rro*(1+alT*(rTemp-25));  
Tgg=gap/2/Rr;

%dTemp: ring temperature increase from id to od (C) in ring fixed reference
%angle_dTemp would be a column containing the angles in degrees and dTemp
%would be a column containing radial ring temperature increase from id to od (in C)
angle_dTemp=linspace(Tgg*180/pi,360-Tgg*180/pi,1000*32+1); %user has to provide this
dTemp=10*ones(1,1000*32+1); %user has to provide this

%=========================FEM input============================%
Nbe=8; %number of elements
Nbnod=Nbe+1; 
Telt=(2*pi-2*Tgg)/(Nbe); 
Npe=120; %number of points within a element for contact calculation
Nr=Nbe*Npe;
Nge=110; %number of contact points within an element for gas force calculation
Ng=Nbe*Nge;

Le=Rr*Telt;
Tnod=(Tgg:(Telt):2*pi-Tgg)+Tg;
Vnod=[Rr*cos(Tnod),Rr*sin(Tnod)];

id=(1:Ng+1)';
shiftstepg=round(Tg/2/pi*Ng);
idg2r=circshift(id,-shiftstepg);
idr2g=circshift(id,shiftstepg);

mr=rhor*2*pi*Rro*Ac;

%========================bore input==============================%
if IsBDist 
    Ampdata=1;  %1 to provide the amplitudes and phases directly or 0 to provide dr distribution
    if Ampdata
        Ab=[91.65 18 11.93 5.07]*1e-6;    %magnitude(m): 0th order, 2nd order, 3rd order, 4th order, etc
        Phib=[0.5642649 1.831025 1.1599458];  %phase(rad): 2nd order, 3rd order, 4th order, etc
    else
        Maxorder=4; %at least 1
        
        thdr=linspace(0,360,1e4)'; %ask for file containing first column of angles in degrees and second 
        dr=zeros(1,1e4)'; %column containing dr in microns
        dr=dr'*1e-6;
        thdr=thdr';
        
        Ab=zeros(1,Maxorder);
            
        if Maxorder>1
            Phib=zeros(1,Maxorder-1);
        else
            Phib=0;
        end
            
        dtheta=2*pi/1000;
        thetaint=linspace(0,2*pi,1001);
        dr=interp1(thdr*pi/180,dr,thetaint);
        fint=dr.*cos(pi/180*thetaint);
        ak=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
        fint=dr.*sin(pi/180*thetaint);
        bk=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
        Ab(1)=sqrt(ak^2+bk^2);
            
        if Maxorder>1
                
           for i=2:Maxorder
                   
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
    beta_p=0.001; %Piston tilt angle, rad
else 
    beta_p=0;
end

off=0; %Piston offset, m

uref=5e-6; 
href=beta_p*Db/2;
if href==0
    href=20e-6;
end

%======================== Oil Data ========================%
zk=0.076642; %zk coefficient for oil viscosity calculation
temp1=1009.863; %temp1 coefficient for oil viscosity calculation, C
temp2=109.047; %temp2 coefficient for oil viscosity calculation, C
rho_oil=850; %Oil density, kg/m^3
hlratio=0.74; %hlratio coefficient for oil viscosity calculation
bta1=2.5; %Beta1 coefficient for oil viscosity calculation
bta2=0.0213; %Beta2 coefficient for oil viscosity calculation
zm0=1;

angle_thetag=linspace(Tgg*180/pi,360-Tgg*180/pi,Ng+1);
angle_theta=linspace(Tgg*180/pi,360-Tgg*180/pi,Nr+1);
anglesg=linspace(0,360,Ng+1);
anglesr=linspace(0,360,Nr+1);
shiftstepr=round(Tg/2/pi*Nr);
shiftstepg=round(Tg/2/pi*Ng);

angle_ringg=linspace(Tgg*180/pi,360-Tgg*180/pi,1000*Nbe+1);

dTemp=interp1(angle_dTemp,dTemp,angle_ringg,'linear','extrap');

%mu_oil: from Thrust side to thrust side in bore fixed reference
%angle_mu_oil would be a column containing the angles in degrees and mu_oil
%would be a column containing oil viscosity in the groove in Pa.s
angle_mu_oil=linspace(0,360,Ng+1); %user has to provide this
mu_oil=1e-6*ones(1,Ng+1); %user has to provide this

mu_oil=interp1(angle_mu_oil,mu_oil,anglesg,'linear','extrap');
mu_oil=circshift(mu_oil,[0 -shiftstepg]);
mu_oil=interp1(anglesg,mu_oil,angle_thetag,'linear','extrap');

%======================liner input================%
El=152.3e9; %liner Young modulus, Pa
nul=0.3; %liner poisson ratio
PR=1; %plateau ratio
sigmap=0.3e-6;  %liner roughness, m

%======================piston input================%
Drldu=94.555e-3;             %upper land reference diameter, m
Drg=86.84e-3;               %ring groove root diameter, m
Drldl=94.555e-3;             %lower land reference diameter, m
hgi=2.05e-3;                %ring groove inner axial height, m
thgt=0;                    %ring groove upper flank angle, rad
thgb=0;                    %ring groove lower flank angle, rad
sigmagt=0.3e-6;            %ring groove upper flank roughness, m
sigmagb=0.3e-6;            %ring groove lower flank roughness, m
nug=0.3;                   %groove poisson ratio
Eg=1.2e11;                 %groove Young modulus, Pa

exp_land_u=210.4e-6;       %radius increase at the upper land during engine operation, m
exp_land_l=177.4e-6;       %radius increase at the lower land during engine operation, m
gr_exp=169.1e-6;           %top ring groove radius increase during engine operation, m

%======================contact input================%
z=6.804;
K=1.198e-4;
A=4.4068e-5;
omega=4;
fc_dry=0.13; %Friction coefficient for the asperity ring/liner and ring/lower plate contact interaction
cfct=1;
Pk_rl=cfct*K*A*2./((1-nul^2)/El+(1-nur^2)/Er);
Pk_rg=cfct*K*A*2./((1-nug^2)/Eg+(1-nur^2)/Er);

y1u=-aui;
y2u=arm-max((arm-auo),(Db-Drldu)/2+Ab(1)-exp_land_u);
lfu=y2u-y1u;
y1l=-ali;
y2l=arm-max((arm-alo),(Db-Drldl)/2+Ab(1)-exp_land_l);
lfl=y2l-y1l;

%======================groove input================%
if IsGDef 
    Ampdatau=1; %1 to provide the amplitudes and phases directly for the upper groove deformation or 0 to provide the axial coordinate distribution
    if Ampdatau
        Agu=[20,0,0]*1e-6*0;   %2nd ring groove upper flank thermal distorsion magnitude: 2nd, 3rd, 4th order, etc
        Phigu=[0,0,0]; %phase(rad): 2nd order, 3rd order, 4th order, etc
    else
        Maxorderu=4; % Number of orders starting from the 2nd one in the Fourier series for the upper groove deformation

        thdrgu=linspace(0,360,1e4)'; %ask for file containing first column of angles in degrees and second 
        dzgu=zeros(1,1e4)'; %column containing dz in microns
        
        dzgu=dzgu'*1e-6;
        thdrgu=thdrgu';
        
        Agu=zeros(1,Maxorderu-1);
        Phigu=zeros(1,Maxorderu-1);
        
        dtheta=2*pi/1000;
        thetaint=linspace(0,2*pi,1001);
        dzgu=interp1(thdrgu*pi/180,dzgu,thetaint,'linear','extrap');
       
        for i=2:Maxorderu
            
            fint=dzgu.*cos(i*pi/180*thetaint);
            ak=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
            fint=dzgu.*sin(i*pi/180*thetaint);
            bk=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
            Agu(i)=sqrt(ak^2+bk^2);
            Phigu(i)=atan2(ak,bk);
        
        end
                
    end
    Ampdatal=1; %1 to provide the amplitudes and phases directly for the lower groove deformation or 0 to provide the axial coordinate distribution
    if Ampdatal
        Agl=[20,0,0]*1e-6*0;   %2nd ring groove lower flank thermal distorsion magnitude; 2nd, 3rd, 4th order, etc
        Phigl=[0,0,0]; %phase(rad): 2nd order, 3rd order, 4th order, etc
    else
        
        Maxorderl=4; % Number of orders starting from the 2nd one in the Fourier series for the lower groove deformation

        thdrgl=linspace(0,360,1e4)'; %ask for file containing first column of angles in degrees and second 
        dzgl=zeros(1,1e4)'; %column containing dz in microns
        
        dzgl=dzgl'*1e-6;
        thdrgl=thdrgl';
        
        Agl=zeros(1,Maxorderl-1);
        Phigl=zeros(1,Maxorderl-1);
        
        dtheta=2*pi/1000;
        thetaint=linspace(0,2*pi,1001);
        dzgl=interp1(thdrgl*pi/180,dzgl,thetaint,'linear','extrap');
        
        for i=2:Maxorderl
            
            fint=dzgl.*cos(i*pi/180*thetaint);
            ak=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
            fint=dzgl.*sin(i*pi/180*thetaint);
            bk=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
            Agl(i)=sqrt(ak^2+bk^2);
            Phigl(i)=atan2(ak,bk);
        
        end
                 
    end
        
else 
    Agu=[0,0,0];   
    Phigu=[0,0,0];
    Agl=[0,0,0];   
    Phigl=[0,0,0];

end

tilt_thu=0; %Groove thermal tilting - upper flank (rad)
tilt_thl=0; %Groove thermal tilting - lower flank (rad)

hrc=hui+thrt*aui+hli+thrb*ali;
hgc=hgi+(Db/2-Drg/2-gr_exp-arm)*(thgt+thgb);
gcl=hgc-hrc;

%=====================lubrication&dynamic input======================%
minoil=sigmap/2; %Minimum oil film thickness to consider for liner lubrication (suggested value: half the surface roughness standard deviation)
oilthreshold=10; %Ratio between oil film thickness on liner by liner surface roughness standard deviation to consider as threshold between parially and fully flooded cases (suggested value: 10)
viscosityfactor=100; %Factor by which the oil viscosity should be divided in regions where we have fuel spots (suggested value: 100)

%Gas pressure in upper region: from Thrust side to thrust side in bore
%fixed reference 
%angle_Pu would be a column containing the angles in degrees and Pu
%would be a column containing the pressure in Bar
angle_Pu=linspace(0,360,Ng+1); %user has to provide this
Pu=0.962*ones(length(angle_Pu),1); %user has to provide this
%Gas pressure in inner region: from Thrust side to thrust side in bore
%fixed reference 
%angle_Pi would be a column containing the angles in degrees and Pi
%would be a column containing the pressure in Bar
angle_Pi=linspace(0,360,Ng+1); %user has to provide this
Pi=0.962*ones(length(angle_Pi),1); %user has to provide this
%Gas pressure in lower region: from Thrust side to thrust side in bore
%fixed reference 
%angle_Pd would be a column containing the angles in degrees and Pd
%would be a column containing the pressure in Bar
angle_Pd=linspace(0,360,Ng+1); %user has to provide this
Pd=0.9067*ones(length(angle_Pd),1); %user has to provide this

Pu=1e5*Pu;
Pi=1e5*Pi;
Pd=1e5*Pd;

Pi=interp1(angle_Pi,Pi,anglesg,'linear','extrap');
Pd=interp1(angle_Pd,Pd,anglesg,'linear','extrap');
Pu=interp1(angle_Pu,Pu,anglesg,'linear','extrap');

Pi=circshift(Pi,[0 -shiftstepg]);
Pi=interp1(anglesg,Pi,angle_thetag,'linear','extrap');
Pd=circshift(Pd,[0 -shiftstepg]);
Pd=interp1(anglesg,Pd,angle_thetag,'linear','extrap');
Pu=circshift(Pu,[0 -shiftstepg]);
Pu=interp1(anglesg,Pu,angle_thetag,'linear','extrap');

Vp=10;  %Piston speed, m/s  
Ap=-0; %Piston acceleration, m^2/s
Fi=-mr*Ap;

%oil film thickness left on liner into the ring: from Thrust side to thrust side in bore
%fixed reference 
%angle_hrls would be a column containing the angles in degrees and hrls
%would be a column containing the oil film thickness in microns
angle_hrls=linspace(0,360,Nr+1); %user has to provide this
hrls=2.5*sigmap*1e6*ones(1,Nr+1);  %user has to provide this

hrls=1e-6*hrls;
hrls=interp1(angle_hrls,hrls,anglesr,'linear','extrap');
hrls=circshift(hrls,[0 -shiftstepr]);
hrls=interp1(anglesr,hrls,angle_theta,'linear','extrap');

%oil film thickness on the upper groove flank: from Thrust side to thrust
%side in bore fixed reference 
%angle_hot would be a column containing the angles in degrees and hot
%would be a column containing the oil film thickness in microns
angle_hot=linspace(0,360,Ng+1); %user has to provide this
hot=0*ones(1,Ng+1); %user has to provide this

hot=1e-6*hot;
hot=interp1(angle_hot,hot,anglesg,'linear','extrap');
hot=circshift(hot,[0 -shiftstepg]);
hot=interp1(anglesg,hot,angle_thetag,'linear','extrap');

%oil film thickness on the lower groove flank: from Thrust side to thrust
%side in bore fixed reference 
%angle_hob would be a column containing the angles in degrees and hob
%would be a column containing the oil film thickness in microns
angle_hob=linspace(0,360,Ng+1); %user has to provide this
hob=0*exp(-(angle_hob-270).^2/2/2/360); %user has to provide this

hob=1e-6*hob;
hob=interp1(angle_hob,hob,anglesg,'linear','extrap');
hob=circshift(hob,[0 -shiftstepg]);
hob=interp1(anglesg,hob,angle_thetag,'linear','extrap'); 

%liner temperature: from Thrust side to thrust side in bore fixed reference
%this temperature distribution will be used to compute oil viscosity
%angle_templ would be a column containing the angles in degrees and templ
%would be a column containing liner temperature in celsius
angle_templ=linspace(0,360,Nr+1); %user has to provide this
templ=150*ones(1,Nr+1);  %user has to provide this

templ=interp1(angle_templ,templ,anglesr,'linear','extrap');
templ=circshift(templ,[0 -shiftstepr]);
templ=interp1(anglesr,templ,angle_theta,'linear','extrap');

%fuel spot: from Thrust side to thrust side in bore fixed reference
%angle_fuel would be a column containing the angles in degrees and isfuel
%would be a column containing fuel spots (1 for existence of fuel)
angle_fuel=linspace(0,360,Nr+1); %user has to provide this
isfuel=zeros(1,Nr+1);  %user has to provide this

isfuel=interp1(angle_fuel,isfuel,anglesr,'linear','extrap');
isfuel=circshift(isfuel,[0 -shiftstepr]);
isfuel=interp1(anglesr,isfuel,angle_theta,'linear','extrap');
isfuel=logical(isfuel);

%=====================fully flooded correlation======================%
b=2.1;
c=0.102;
d=-0.04;
e=0.90229;

%=====================deterministic correlation======================%
Ph=3.9856e7; 
KOCR=2.5230; 
muUot=1.5e-2;

ap=0.039342;
cp1=-2.965331;
cp2=1.499148;
F0=0.098577;
cf1=-0.383954;
cf2=0.138443;
rw=0.7e-3; %Ring width to use for the determinstic hydrodynamic correlation, m

%=====================prepare input matrix======================%
inputs.pistonSM.beta_p=beta_p;
inputs.pistonSM.off=off;

inputs.bore.Db=Db;
inputs.bore.Ab=Ab;
inputs.bore.Phib=Phib;
inputs.bore.templ=templ;

inputs.oil.zk=zk;
inputs.oil.temp1=temp1;
inputs.oil.temp2=temp2;
inputs.oil.rho_oil=rho_oil;
inputs.oil.hlratio=hlratio;
inputs.oil.bta1=bta1;
inputs.oil.bta2=bta2;
inputs.oil.zm0=zm0;
inputs.oil.mu_oil=mu_oil;

inputs.ring.auo=auo;
inputs.ring.aui=aui;
inputs.ring.alo=alo;
inputs.ring.ali=ali;
inputs.ring.huo=huo;
inputs.ring.hui=hui;
inputs.ring.hlo=hlo;
inputs.ring.hli=hli;
inputs.ring.thrt=thrt;           
inputs.ring.thrb=thrb;           
inputs.ring.rb1=rb1;
inputs.ring.rb2=rb2;
inputs.ring.rbn=rbn;
inputs.ring.a10=a10;
inputs.ring.a11=a11;
inputs.ring.a12=a12;
inputs.ring.a20=a20;
inputs.ring.a21=a21;
inputs.ring.a22=a22;
inputs.ring.arm=arm;
inputs.ring.gap=gap;
inputs.ring.rhor=rhor;
inputs.ring.Er=Er;
inputs.ring.alT=alT;
inputs.ring.nur=nur;
inputs.ring.Gr=Gr;
inputs.ring.Ac=Ac;
inputs.ring.Izz=Izz;
inputs.ring.Iyy=Iyy;
inputs.ring.Izr=Izr;
inputs.ring.Ip=Ip;
inputs.ring.Jt=Jt;
inputs.ring.alp=alp;
inputs.ring.Rr=Rr;
inputs.ring.Le=Le;
inputs.ring.mr=mr;
inputs.ring.Tg=Tg;
inputs.ring.rTemp=rTemp;
inputs.ring.dTemp=dTemp;
inputs.ring.idr2g=idr2g;
inputs.ring.idg2r=idg2r;

inputs.piston.Drldu=Drldu;
inputs.piston.Drg=Drg;
inputs.piston.Drldl=Drldl;
inputs.piston.hgi=hgi;
inputs.piston.thgt=thgt;
inputs.piston.thgb=thgb;
inputs.piston.sigmag1t=sigmagt;
inputs.piston.sigmag1b=sigmagb;
inputs.piston.nug=nug;
inputs.piston.Eg=Eg;
inputs.piston.exp_land_u=exp_land_u;
inputs.piston.exp_land_l=exp_land_l;
inputs.piston.gr_exp=gr_exp;

inputs.liner.El=El;
inputs.liner.nul=nul;
inputs.liner.sigmap=sigmap;
inputs.liner.PR=PR;

inputs.FEM.Nbe=Nbe;
inputs.FEM.Nbnod=Nbnod;
inputs.FEM.Telt=Telt;
inputs.FEM.Npe=Npe;
inputs.FEM.Nr=Nr;
inputs.FEM.Le=Le;
inputs.FEM.Tnod=Tnod;
inputs.FEM.Vnod=Vnod;
inputs.FEM.uref=uref;
inputs.FEM.href=href;
inputs.FEM.Nge=Nge;
inputs.FEM.Ng=Ng;

inputs.contact.z=z;
inputs.contact.Pk_rl=Pk_rl;
inputs.contact.Pk_rg=Pk_rg;
inputs.contact.omega=omega;
inputs.contact.fc_dry=fc_dry;
inputs.contact.cfct=cfct;

inputs.groove.y1u=y1u;
inputs.groove.y2u=y2u;
inputs.groove.lfu=lfu;
inputs.groove.y1l=y1l;
inputs.groove.y2l=y2l;
inputs.groove.lfl=lfl;
inputs.groove.gcl=gcl;
inputs.groove.Agu=Agu;
inputs.groove.Phigu=Phigu;
inputs.groove.Agl=Agl;
inputs.groove.Phigl=Phigl;
inputs.groove.tilt_thu=tilt_thu;
inputs.groove.tilt_thl=tilt_thl;

inputs.engine.Vp=Vp;
inputs.engine.Ap=Ap;
inputs.ring.Fi=Fi;
inputs.engine.Pu=Pu;
inputs.engine.Pd=Pd;
inputs.engine.Pi=Pi;

inputs.liner.PR=PR;
inputs.liner.Ph=Ph;
inputs.liner.KOCR=KOCR;
inputs.liner.ap=ap;
inputs.liner.cp1=cp1;
inputs.liner.cp2=cp2;
inputs.liner.F0=F0;
inputs.liner.cf1=cf1;
inputs.liner.cf2=cf2;
inputs.liner.rw=rw;
inputs.liner.muUot=muUot;
inputs.liner.b=b;
inputs.liner.c=c;
inputs.liner.d=d;
inputs.liner.e=e;
inputs.lubrication.hrls=hrls;
inputs.lubrication.isfuel=isfuel;
inputs.lubrication.hot=hot;
inputs.lubrication.hob=hob;
inputs.lubrication.minoil=minoil;
inputs.lubrication.oilthreshold=oilthreshold;
inputs.lubrication.viscosityfactor=viscosityfactor;

%=============freeshape================%
if IsFreeshape
    IsFreeshaped=0; %1 to provide the free shape coordinates and 0 to provide the free shape curvature
    if IsFreeshaped
        %ask for free shape data with first column containing angles in degrees and
        %the second one containing the freeshape in millimeters
        FS_raw=dlmread('ringfreeshape.txt'); %user has to provide this
        Tri=FS_raw(:,1)*pi/180;
        rfsim=FS_raw(:,2)*1e-3;

        Inod=1:1000:1000*Nbe+1;
        dTrfs=(Tri(end)-Tri(1))/(1000*Nbe);
        Trfs=(Tri(1):dTrfs:Tri(length(Tri)))';
        Teltfs=(Trfs(end)-Trfs(1))/Nbe;
        Tnodfs=Trfs(Inod);
    
        rfsi=interp1(Tri,rfsim,Trfs);
        [Ccs]=Cs_3rd_cont_n4g(Tnodfs);
        [rnodfs]=ls_rfsi_n4g(rfsi,Nbnod,Teltfs,Ccs);
        [Kapfs,~]=Kfs_cal_n4g(Trfs,rfsi,rnodfs,Nbe,1000,Teltfs);
    else
    
        %ask for free shape curvature with first column containing angles in degrees and
        %the second one containing the freeshape curvature in 1/m
    
        Kapfs=dlmread('C:\Users\YangLiu\Dropbox (MIT)\Aziz Perso\papers mit\RA\Thesis\Texts\Kapfs_from_theo_pr_TC.txt');
        Tri=Kapfs(:,1)*pi/180;
        Kapfs=interp1(Tri,Kapfs(:,2),linspace(Tgg,2*pi-Tgg,1000*Nbe+1)','linear','extrap');
    end
else
    % ask for text file containing two columns
    %first column containing angles in degrees
    %second column containing pressure times ring height (axial width) i.e.
    %2.Lt/Dn*NormalizedPr
    PressInput=dlmread('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/theo_pr_TC1.txt'); %user has to provide this
    
    Tpo=PressInput(:,1);
    Pmo=PressInput(:,2);
    [~,~,Kapfs]=fscal_n4g(Tpo,Pmo,Rr,Er,Izr,gap,Nbe,1000);
end

hr=hui+thrt*aui+hli+thrb*ali;
ar=arm+max(ali,aui);

%=====================solver======================%
Newton_tol=1e-6; %Tolerance for Newton-Raphson algorithm convergence (suggested value: 1e-6)
ItMax=100; %Maximum number of Newton-Raphson algorithm iterations (suggested value: 100)

Ug1=zeros(8*Nbnod,1);
Ug1(7:8:end-1)=alp;

F1=Finitial_fs_vg(Ug1,Er,Gr,Iyy,Izz,Jt,alp,Rr,Telt,Nbe,1000,Kapfs,dTemp,alT,hr,ar);

inputs.ring.F1=F1;

[Kg_nd]=Kg_ring_curved_nd_vg(Er,Gr,Iyy,Jt,Rr,Telt,Nbe);

inputs.FEM.Kg_nd=Kg_nd;

UgIC_nd=zeros(8*Nbnod,1);
[ybn]=bore_dist_static_vg(inputs);
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
[zgnu,zgnl]=groove_dist_static_vg(inputs);

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
[NLS_0,yr,zr,alr,hmin,fl,z0,~,~,~]=Static_Error_vg(Ug_nd,F1,inputs);
NLS_k=NLS_0;
j=0;
while j<=ItMax&&(judge1||judge2)
    j=j+1;
    [Ja_stac]=Jacobian_Estac_vg(Ug_nd,inputs);
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
            [NLS_k,yr,zr,alr,hmin,fl,z0,~,~,~]=Static_Error_vg(Ug_nd,F1,inputs);
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
            [NLS_k,yr,zr,alr,hmin,fl,z0,~,~,~]=Static_Error_vg(Ug_nd,F1,inputs);
            Err=norm(NLS_k); 
            dErr=Erro-Err;
            itGB=itGB+1;
        end
    end
   
    judge1=norm(Nstep)>Newton_tol;
    judge2=norm(NLS_k)>Newton_tol;
    
    norm_NLS_k(j)=norm(NLS_k);
    norm_step(j)=norm(Nstep);
    
    figure(1)
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
  
 [ybn]=bore_dist_static_vg(inputs);
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

 [zgnu,zgnl]=groove_dist_static_vg(inputs);
 
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
 yring=zeros(1,Nr+1);
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
     yring(ind)=yrl;
     
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
 
  display(['Ft(CBM) = ', num2str(sum(fl*Rr/Nr)),'N']) %tengential load
 theer=linspace(Tgg*180/pi,360-Tgg*180/pi,Nr+1)';
 hmin=hmin';
 hmin=1e6*max(0,hmin);
 z0=1e6*z0';
 zr=1e6*zr';
 yr=1e6*yr';
 alr=180/pi*alr';
 uoc=1e6*uoc';
 uoc=max(0,uoc);
 uic=1e6*uic';
 uic=max(0,uic);
 loc=1e6*loc';
 loc=max(0,loc);
 lic=1e6*lic';
 lic=max(0,lic);
 
  figure(1)
 clf
 set(gcf,'Position',[314    58   703   916])
 subplot(211)
 hold on
 plot(theer(1:Npe:end),hmin(1:Npe:end),'ro')
 plot(theer,hmin,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Ring Liner Clearance (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 subplot(212)
 hold on
 plot(theer(1:Npe:end),yr(1:Npe:end),'ro')
 plot(theer,yr,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Radial coordinate % nominal radius (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 
 figure(2)
 subplot(311)
 hold on
 plot(theer(1:Npe:end),zr(1:Npe:end),'ro')
 plot(theer,zr,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Ring Netrual Axis Lift (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 subplot(312)
 hold on
 plot(theer(1:Npe:end),alr(1:Npe:end),'ro')
 plot(theer,alr,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Ring Twist Angle (degree)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 subplot(313)
 hold on
 plot(theer(1:Npe:end),z0(1:Npe:end),'ro')
 plot(theer,z0,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Min clearance pt location % centroid (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 
 figure(3)
 clf
 set(gcf,'Position',[144 64  1302 916])
 subplot(221)
 hold on
 plot(theer(1:Npe:end),uoc(1:Npe:end),'ro')
 plot(theer,uoc,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Upper OD Clearance (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 subplot(222)
 hold on
 plot(theer(1:Npe:end),uic(1:Npe:end),'ro')
 plot(theer,uic,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Upper ID Clearance (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 subplot(223)
 hold on
 plot(theer(1:Npe:end),loc(1:Npe:end),'ro')
 plot(theer,loc,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Lower OD Clearance (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 subplot(224)
 hold on
 plot(theer(1:Npe:end),lic(1:Npe:end),'ro')
 plot(theer,lic,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Lower ID Clearance (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 
 mat_res_1=[theer hmin zr yr alr z0 uoc uic loc lic];
 
 temp=inputs.bore.templ;
 dt=1;
 sigmagt=inputs.piston.sigmag1t;
 sigmagb=inputs.piston.sigmag1b;
 Ph_OCR=inputs.liner.Ph;
 fc=inputs.contact.fc_dry;
 isfuela=inputs.lubrication.isfuel;
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
    [~,~,fhl,ffl,mhl]=ring_liner_hydro_vg(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yre,alre,hrle,KOCR,Kpe,Kfe,Ph_OCR,ap,F0,b,c,d,e,rw,muUot,templ,isfuele,Vp,inputs.oil,sigmap,oilthreshold,viscosityfactor);
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
 Nplot=min(Ng,Nr);
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
 fy=-fy;
 fl=fclglob+fhlglob;
fl=-fl;

[ybn]=bore_dist_static_vg(inputs);
 for j=1:Nbe
     ind=1+(j-1)*Npe:1+j*Npe;
     id=3*(j-1)+1:3*(j+1);
     yb=(ybn(id))'*N*uref;
     
     ybr(ind)=yb;
 end
 
figure(4)
clf
set(gcf,'Position',[144 64  1302 916])
subplot(221)
hold on
plot(thee(1:min(Nge,Npe):end),fy(1:min(Nge,Npe):end),'ro')
plot(thee,fy,'m','linewidth',1.5)
xlabel('Circumferential Direction (degree)','FontSize',12)
ylabel('Radial Force (N/m)','FontSize',12)
set(gca,'XTick',[0:90:360])
xlim([0 360])
hold off
grid
set(gca,'FontSize',16);
subplot(222)
hold on
plot(thee(1:min(Nge,Npe):end),fl(1:min(Nge,Npe):end),'ro')
plot(thee,fl,'m','linewidth',1.5)
xlabel('Circumferential Direction (degree)','FontSize',12)
ylabel('RIng Liner Radial Force (N/m)','FontSize',12)
set(gca,'XTick',[0:90:360])
xlim([0 360])
hold off
grid
set(gca,'FontSize',16);
subplot(223)
hold on
plot(thee(1:min(Nge,Npe):end),fz(1:min(Nge,Npe):end),'ro')
plot(thee,fz,'m','linewidth',1.5)
xlabel('Circumferential Direction (degree)','FontSize',12)
ylabel('Axial Force (N/m)','FontSize',12)
set(gca,'XTick',[0:90:360])
xlim([0 360])
hold off
grid
set(gca,'FontSize',16);
subplot(224)
hold on
plot(thee(1:min(Nge,Npe):end),m(1:min(Nge,Npe):end),'ro')
plot(thee,m,'m','linewidth',1.5)
xlabel('Circumferential Direction (degree)','FontSize',12)
ylabel('Twist moment (N)','FontSize',12)
set(gca,'XTick',[0:90:360])
xlim([0 360])
hold off
grid
set(gca,'FontSize',16);

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

dx=1/min(Nge,Npe);
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

for i=1:Nbe-1
    
    Uey=Uov(3*i-2:3*i+3);
    Uez=Uezn(3*i-2:3*i+3);
    Uea=Uean(3*i-2:3*i+1);
    
    ov(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uey')*N(:,1:end-1);
    dov(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uey')*dN(:,1:end-1)/Telt;
    d2ov(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uey')*d2N(:,1:end-1)/Telt^2;
    
    dz(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uez')*dN(:,1:end-1)/Telt;
    d2z(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uez')*d2N(:,1:end-1)/Telt^2;
    
    alpha(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uea')*Nt(:,1:end-1);
    dalpha(1+(i-1)*min(Nge,Npe):i*min(Nge,Npe))=(Uea')*dNt(:,1:end-1)/Telt;
    
end

Uey=Uov(end-5:end);
ov(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uey')*N;
dov(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uey')*dN/Telt;
d2ov(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uey')*d2N/Telt^2;

Uez=Uezn(end-5:end);
dz(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uez')*dN/Telt;
d2z(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uez')*d2N/Telt^2;

Uea=Uean(end-3:end);
alpha(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uea')*Nt;
dalpha(1+(Nbe-1)*min(Nge,Npe):1+Nbe*min(Nge,Npe))=(Uea')*dNt/Telt;

kzz=(1/Rr-(ov+d2ov)/Rr^2).*sin(alpha+d2z/Rr);
kzz=kzz';

Kapfs=interp1(linspace(Tgg,2*pi-Tgg,1000*Nbe+1)',Kapfs,linspace(Tgg,2*pi-Tgg,Nbe*min(Nge,Npe)+1)','linear','extrap');

kapfsz=Kapfs*sin(alp);

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

Pm=fy;
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
 
 figure(7)
set(gcf,'Position',[144 64  1302 916])
hold on
plot(thee(1:min(Nge,Npe):end),Mfinal(1:min(Nge,Npe):end),'ro')
plot(thee,Mfinal,'m','linewidth',1.5)
xlabel('Circumferential Direction (degree)','FontSize',12)
ylabel('Bending moment (Nm)','FontSize',12)
set(gca,'XTick',[0:90:360])
xlim([0 360])
hold off
grid
 set(gca,'FontSize',16);

 figure(8)
set(gcf,'Position',[144 64  1302 916])
hold on
plot(thee,S_uo,'b','linewidth',1.5)
plot(thee,S_ui,'k','linewidth',1.5)
plot(thee,S_lo,'--','Color',[1 0 0],'linewidth',1.5)
plot(thee,S_li,'--','Color',[1 1 0],'linewidth',1.5)
xlabel('Circumferential Direction (degree)','FontSize',12)
ylabel('Stress (MPa)','FontSize',12)
legend('Stress at Upper OD','Stress at Upper ID','Stress at Lower OD','Stress at Lower ID','Location','Best')
set(gca,'XTick',[0:90:360])
xlim([0 360])
hold off
grid
 set(gca,'FontSize',16);

 mat_res_2=[thee fy fl fz m Mfinal S_uo S_ui S_lo S_li];
 
 hminv=1e-6*hmin;
dtg=180*gap/pi/2/Rr+0.2;

if IsBDist
theeplotv=linspace(0,360,Nplot+1);
theev=linspace(0,360,Nr+1)';
 
rbr=mag1*interp1(theev,ybr',theeplotv)+Db/2;
format long
rbmax=max(abs(interp1(theev,ybr',theeplotv))); 
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
rbr=mag1*interp1(theev,ybr',theeplotv)+Db/2;
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
fyp=fy;
flp=fl;

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
flp=rmax/2/flmax*flp.*(abs(flp)<fstop)+rmax*sign(flp).*(abs(flp)>=fstop);

numplot=1000;

figure(5)
clf
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

figure(6)
clf
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
fyp=fy;
flp=fl;

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
flp=rmax/2/flmax*flp.*(abs(flp)<fstop)+rmax*sign(flp).*(abs(flp)>=fstop);

numplot=1000;

figure(5)
clf
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
figure(6)
clf
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
