function rdt_conformability_video_test_g()

close all

%this code performs the calculation for the conformability study for different ring gap positions equally spaced between 0 and 359 degrees

addpath(genpath(pwd))

%=============features================================%
IsBDist=0; %1 to take into account bore thermal distortion, 0 otherwise
IsPTilt=0; %1 to take into account piston dynamic tilt, 0 otherwise
IsGDef=0; %1 to take into account groove thermal distortion, 0 otherwise
IsFreeshape=0;  %1 to provide free shape or its curvature and 0 to provide radial pressure distribution for the circular shape
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
gap=0.3314e-3; %ring gap size when closed inside bore, m

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
id=(1:Ng+1)';
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

angle_ringg=linspace(Tgg*180/pi,360-Tgg*180/pi,1000*Nbe+1);

dTemp=interp1(angle_dTemp,dTemp,angle_ringg,'linear','extrap');

%mu_oil: from Thrust side to thrust side in bore fixed reference
%angle_mu_oil would be a column containing the angles in degrees and mu_oil
%would be a column containing oil viscosity in the groove in Pa.s
angle_mu_oil=linspace(0,360,Ng+1); %user has to provide this
mu_oilglob=1e-6*ones(1,Ng+1); %user has to provide this

mu_oilglob=interp1(angle_mu_oil,mu_oilglob,anglesg,'linear','extrap');

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
%angle_Pu would be a column containing the angles in degrees and Puglob
%would be a column containing the pressure in Bar
angle_Pu=linspace(0,360,Ng+1); %user has to provide this
Puglob=0.962*ones(length(angle_Pu),1); %user has to provide this
%Gas pressure in inner region: from Thrust side to thrust side in bore
%fixed reference 
%angle_Pi would be a column containing the angles in degrees and Piglob
%would be a column containing the pressure in Bar
angle_Pi=linspace(0,360,Ng+1); %user has to provide this
Piglob=0.962*ones(length(angle_Pi),1); %user has to provide this
%Gas pressure in lower region: from Thrust side to thrust side in bore
%fixed reference 
%angle_Pd would be a column containing the angles in degrees and Pdglob
%would be a column containing the pressure in Bar
angle_Pd=linspace(0,360,Ng+1); %user has to provide this
Pdglob=0.9067*ones(length(angle_Pd),1); %user has to provide this

Puglob=1e5*Puglob;
Piglob=1e5*Piglob;
Pdglob=1e5*Pdglob;

Piglob=interp1(angle_Pi,Piglob,anglesg,'linear','extrap');
Pdglob=interp1(angle_Pd,Pdglob,anglesg,'linear','extrap');
Puglob=interp1(angle_Pu,Puglob,anglesg,'linear','extrap');

Vp=10;  %Piston speed, m/s  
Ap=-0; %Piston acceleration, m^2/s
Fi=-mr*Ap;

%oil film thickness left on liner into the ring: from Thrust side to thrust side in bore
%fixed reference 
%angle_hrls would be a column containing the angles in degrees and hrlsglob
%would be a column containing the oil film thickness in microns
angle_hrls=linspace(0,360,Nr+1); %user has to provide this
hrlsglob=2.5*sigmap*1e6*ones(1,Nr+1);  %user has to provide this

hrlsglob=1e-6*hrlsglob;
hrlsglob=interp1(angle_hrls,hrlsglob,anglesr,'linear','extrap');

%oil film thickness on the upper groove flank: from Thrust side to thrust
%side in bore fixed reference 
%angle_hot would be a column containing the angles in degrees and hotglob
%would be a column containing the oil film thickness in microns
angle_hot=linspace(0,360,Ng+1); %user has to provide this
hotglob=0*ones(1,Ng+1); %user has to provide this

hotglob=1e-6*hotglob;
hotglob=interp1(angle_hot,hotglob,anglesg,'linear','extrap');

%oil film thickness on the lower groove flank: from Thrust side to thrust
%side in bore fixed reference 
%angle_hob would be a column containing the angles in degrees and hobglob
%would be a column containing the oil film thickness in microns
angle_hob=linspace(0,360,Ng+1); %user has to provide this
hobglob=0*exp(-(angle_hob-270).^2/2/2/360); %user has to provide this

hobglob=1e-6*hobglob;
hobglob=interp1(angle_hob,hobglob,anglesg,'linear','extrap');

%liner temperature: from Thrust side to thrust side in bore fixed reference
%this temperature distribution will be used to compute oil viscosity
%angle_templ would be a column containing the angles in degrees and templglob
%would be a column containing liner temperature in celsius
angle_templ=linspace(0,360,Nr+1); %user has to provide this
templglob=150*ones(1,Nr+1);  %user has to provide this

templglob=interp1(angle_templ,templglob,anglesr,'linear','extrap');

%fuel spot: from Thrust side to thrust side in bore fixed reference
%angle_fuel would be a column containing the angles in degrees and isfuelglob
%would be a column containing fuel spots (1 for existence of fuel)
angle_fuel=linspace(0,360,Nr+1); %user has to provide this
isfuelglob=zeros(1,Nr+1);  %user has to provide this

isfuelglob=interp1(angle_fuel,isfuelglob,anglesr,'linear','extrap');

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

inputs.oil.zk=zk;
inputs.oil.temp1=temp1;
inputs.oil.temp2=temp2;
inputs.oil.rho_oil=rho_oil;
inputs.oil.hlratio=hlratio;
inputs.oil.bta1=bta1;
inputs.oil.bta2=bta2;
inputs.oil.zm0=zm0;

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

inputs.ring.rTemp=rTemp;
inputs.ring.dTemp=dTemp;

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

%=====================solver======================%;
Newton_tol=1e-6; %Tolerance for Newton-Raphson algorithm convergence (suggested value: 1e-6)
ItMax=100; %Maximum number of Newton-Raphson algorithm iterations (suggested value: 100)
numTg=360; %Number of ring gap positions equally spaced between 0 and 359 to consider

Ug1=zeros(8*Nbnod,1);
Ug1(7:8:end-1)=alp;

F1=Finitial_fs_vg(Ug1,Er,Gr,Iyy,Izz,Jt,alp,Rr,Telt,Nbe,1000,Kapfs,dTemp,alT,hr,ar);

Kapfs=interp1(linspace(Tgg,2*pi-Tgg,1000*Nbe+1)',Kapfs,linspace(Tgg,2*pi-Tgg,Nbe*min(Nge,Npe)+1)','linear','extrap');

kapfsz=Kapfs*sin(alp);

inputs.ring.F1=F1;

[Kg_nd]=Kg_ring_curved_nd_vg(Er,Gr,Iyy,Jt,Rr,Telt,Nbe);

inputs.FEM.Kg_nd=Kg_nd;

UgIC_nd=zeros(8*Nbnod,1);
 
 Nplot=min(Ng,Nr);
 res_Tg=linspace(0,359,numTg)';

 res_theer=zeros(Nr+1,numTg);
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
 
 res_thee=zeros(Nplot+1,numTg);
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
jframe
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

inputs.ring.idr2g=idr2g;
inputs.ring.idg2r=idg2r;
inputs.ring.Tg=Tg;
inputs.FEM.Tnod=Tnod;
inputs.FEM.Vnod=Vnod;

inputs.oil.mu_oil=mu_oil;

inputs.engine.Pu=Pu;
inputs.engine.Pd=Pd;
inputs.engine.Pi=Pi;

inputs.lubrication.hrls=hrls;
inputs.lubrication.hot=hot;
inputs.lubrication.hob=hob;
inputs.lubrication.isfuel=isfuel;

inputs.bore.templ=templ;

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
    if mod(j,10)==0
        j
    end
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
    
end
j

hmin=max(0,hmin);

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
 
  display(['Ft(CBM) = ', num2str(sum(fl*Rr/Nr)),'N']) %tengential load

 theer=linspace(Tgg*180/pi,360-Tgg*180/pi,Nr+1)';

 uoc=max(0,uoc);
 uic=max(0,uic);
 loc=max(0,loc);
 lic=max(0,lic);
 
 res_theer(:,jframe)=theer;
 res_hmin(:,jframe)=hmin'*1e6;
 res_zr(:,jframe)=zr'*1e6;
 res_yr(:,jframe)=yr'*1e6;
 res_alr(:,jframe)=alr'/pi*180;
 res_z0(:,jframe)=z0'*1e6;
 res_uoc(:,jframe)=uoc'*1e6;
 res_uic(:,jframe)=uic'*1e6;
 res_loc(:,jframe)=loc'*1e6;
 res_lic(:,jframe)=lic'*1e6;
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
 
 res_thee(:,jframe)=thee;
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

 S_lo=[S_lo(1:end-1);flipud(S_lo2)];
 rlo=sqrt(alo^2+hlo^2);
 S_lo=S_lo+abs(Er*Iyy*(kzz-kapfsz))*hlo/Iyr.*(fz>=0)-abs(Er*Iyy*(kzz-kapfsz))*hlo/Iyr.*(fz>=0)+rlo*Gr*(dz'/Rr^2-dalpha'/Rr);

 S_ui=[S_ui(1:end-1);flipud(S_ui2)];
 rui=sqrt(aui^2+hui^2);
 S_ui=S_ui-abs(Er*Iyy*(kzz-kapfsz))*hui/Iyr.*(fz>=0)+abs(Er*Iyy*(kzz-kapfsz))*hui/Iyr.*(fz>=0)+rui*Gr*(dz'/Rr^2-dalpha'/Rr);

 S_li=[S_li(1:end-1);flipud(S_li2)];
 rli=sqrt(ali^2+hli^2);
 S_li=S_li+abs(Er*Iyy*(kzz-kapfsz))*hli/Iyr.*(fz>=0)-abs(Er*Iyy*(kzz-kapfsz))*hli/Iyr.*(fz>=0)+rli*Gr*(dz'/Rr^2-dalpha'/Rr);
 
 res_Mfinal(:,jframe)=Mfinal;
 res_S_uo(:,jframe)=S_uo*1e-6;
 res_S_lo(:,jframe)=S_lo*1e-6;
 res_S_ui(:,jframe)=S_ui*1e-6;
 res_S_li(:,jframe)=S_li*1e-6;

 [ybn]=bore_dist_static_vg(inputs);
 for j=1:Nbe
     ind=1+(j-1)*Npe:1+j*Npe;
     id=3*(j-1)+1:3*(j+1);
     yb=(ybn(id))'*N*uref;
     
     ybr(ind)=yb;
 end
 res_ybr(:,jframe)=ybr';
 
end
 
 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_Tg_ng.txt',res_Tg,'delimiter','\t','precision','%12.8f','newline','pc');
 
 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_theer_ng.txt',res_theer,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_hmin_ng.txt',res_hmin,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_zr_ng.txt',res_zr,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_yr_ng.txt',res_yr,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_alr_ng.txt',res_alr,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_z0_ng.txt',res_z0,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_uoc_ng.txt',res_uoc,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_uic_ng.txt',res_uic,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_loc_ng.txt',res_loc,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_lic_ng.txt',res_lic,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_ybr_ng.txt',res_ybr,'delimiter','\t','precision','%12.8f','newline','pc');

 
 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_thee_ng.txt',res_thee,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_fy_ng.txt',res_fy,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_fl_ng.txt',res_fl,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_fz_ng.txt',res_fz,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_m_ng.txt',res_m,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_Mfinal_ng.txt',res_Mfinal,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_S_uo_ng.txt',res_S_uo,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_S_lo_ng.txt',res_S_lo,'delimiter','\t','precision','%12.8f','newline','pc');

 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_S_ui_ng.txt',res_S_ui,'delimiter','\t','precision','%12.8f','newline','pc');
 
 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_S_li_ng.txt',res_S_li,'delimiter','\t','precision','%12.8f','newline','pc');


 dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/util_data_ng.txt',util_data,'delimiter','\t','precision','%12.8f','newline','pc');

end