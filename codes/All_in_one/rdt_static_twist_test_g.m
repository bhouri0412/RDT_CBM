function rdt_static_twist_test_g()

close all

%this code performs the static twist calculations under FixOD or FixID Constraint

% it outputs the node's circumferential location (in microns), ring-liner
% clearance (microns), twist angle (degree), lower OD clearance (microns),
% lower ID clearance (microns)
%radial force (N/m), axial force (N/m) and twist moment (N) 
%all the results are contained in the variable mat_result

addpath(genpath(pwd))

Tg=0/180*pi;

IsBDist=0; 
Ampdata=0;
IsPTilt=0; 
IsGDef=0;
%=============features================================%
IsFreeshape=0; %1 if user provides free shape coordinates, 0 if provides free shape curvature
inputs.features.IsBDist=IsBDist;
inputs.features.IsPTilt=IsPTilt;
inputs.features.IsGDef=IsGDef;
inputs.features.IsFreeshape=IsFreeshape;

%============================ring input=======================%
Db=95.25e-3;  %bore diameter, m
auo= 1.909325626204240e-3;   %upper OD width, m
aui= 2.090674373795760e-3;   %upper ID width, m
alo= 1.909325626204240e-3;   %lower OD width, m
ali= 2.090674373795760e-3;   %lower ID width, m
huo= 1.050443159922929e-3;   %upper OD height, m
hui= 1.050443159922929e-3;  %upper ID height, m
hlo= 1.149556840077072e-3;  %lower OD height, m
hli= 0.949556840077072e-3;  %lower ID height, m
thrt= 0*pi/180;           %ring upper flank angle, rad
thrb= 2.862405226111749*pi/180;           %ring lower flank angle, rad
 
rb1=0.7e-3;    %the lower edge width, m
rb2=0.7e-3;     %the upper edge width, m
rbn= 0.050443159922928e-3;  %running face minimum point axial location, m
a10=0; 
a11=0; %Linear cofficient for lower edge shape factor:
a12=50;        %the lower edge shape factor, 1/m
a20=0;
a21=0; %Linear cofficient for upper edge shape factor:
a22=50;       %the upper edge shape factor, 1/m
arm= 2.109325626204241e-3; %running face minimum point width, m
gap=0.314e-3;  %ring gap size when closed, m
 
Ac= 8.65e-6;   %area of cross-section, m^2
Izz= 12.209622500233360e-12; %principal moment of inertial in plane, m^4
Iyy= 3.149781867139798e-12; %principal moment of inertial out of plane, m^4
Izr= 12.201364399486193e-12; %moment of inertial in plane Iz, m^4
Ip=Izz+Iyy;   
Jt= 8.589953839799971e-12;   %torsional factor, m^4
alp= -1.730088983692877/180*pi;  %principal angle, rad
isi=0*(alp<=0)+1*(alp>0);

rhor=7820;       %ring density, kg/m^3
Er=100e9;        %ring Young's modulus, Pa
nur=0.3;         %ring Poisson Ratio
Gr=Er/(2*(1+nur));
alT=1.1e-5;

rTemp=25;

Rro=Db/2-arm;
Rr=Rro*(1+alT*(rTemp-25));  

%=========================FEM input============================%
Nbe=32; %number of elements
Nbnod=Nbe+1; 
Tgg=gap/2/Rr;
Telt=(2*pi-2*Tgg)/(Nbe);
Npe=1000; %number of points within a element for contact calculation
Nr=Nbe*Npe;
Nge=1;
Ng=Nbe*Nge;

Le=Rr*Telt;
Tnod=(Tgg:(Telt):2*pi-Tgg)+Tg;
Vnod=[Rr*cos(Tnod),Rr*sin(Tnod)];

id=(1:Ng+1)';
shiftstep=round(Tg/2/pi*Ng);
idg2r=circshift(id,-shiftstep);
idr2g=circshift(id,shiftstep);

mr=rhor*2*pi*Rro*Ac;

%========================bore input==============================%
if IsBDist 
    if Ampdata
        Maxorder=4; 
        Ab=[3*91.65 18 11.93 5.07]*1e-6;   
        Phib=[0.5642649 1.831025 1.1599458];  
    else
        
        Maxorder=4; 
        thdr=linspace(0,360,1e4);
        dr=zeros(1,1e4);
        
        Ab=zeros(1,Maxorder);
            
        if Maxorder>1
            Phib=zeros(1,Maxorder-1);
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

if IsPTilt
    beta_p=0.001;
else 
    beta_p=0;
end

off=0;

uref=5e-6;  
href=beta_p*Db/2;
if href==0
    href=20e-6;
end

fref=100;

%======================liner input================%
El=152.3e9; %liner Young modulus, Pa
nul=0.3; %liner poisson ratio
PR=1; %plateau ratio
sigmap=0.3e-6;  %liner roughness, m 

Drldu=Db;           
Drg=Db-2*(arm+aui); 
Drldl=Db;     
hgi=2*(hui+hli);      
thgt=0;  
sigmagt=0.4e-6;          
exp_land_u=0;     
exp_land_l=0;    
gr_exp=0; 
thgb=thrb;
sigmagb=0.4e-6;   
nug=0.3; 
Eg=1.2e11;

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

if IsGDef 
    Ampdatau=1;
    if Ampdatau
        Maxorderu=4; 
        Agu=[20,0,0]*1e-6*0; 
        Phigu=[0,0,0];
        
    else
        
        Maxorderu=4;
        
        thdrgu=linspace(0,360,1e4);
        drgu=zeros(1,1e4);
        
        Agu=zeros(1,Maxorderu-1);
        Phigu=zeros(1,Maxorderu-1);
        
        dtheta=2*pi/1000;
        thetaint=linspace(0,2*pi,1001);
        drgu=interp1(thdrgu*pi/180,drgu,thetaint,'linear','extrap');
       
        for i=2:Maxorderu
            
            fint=drgu.*cos(i*pi/180*thetaint);
            ak=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
            fint=drgu.*sin(i*pi/180*thetaint);
            bk=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
            Agu(i)=sqrt(ak^2+bk^2);
            Phigu(i)=atan2(ak,bk);
        
        end
                
    end
    Ampdatal=1;
    if Ampdatal
        Maxorderl=4;
        Agl=[20,0,0]*1e-6*0; 
        Phigl=[0,0,0];

    else
        
        Maxorderl=4;

        thdrgl=linspace(0,360,1e4);
        drgl=zeros(1,1e4);
        
        Agl=zeros(1,Maxorderl-1);
        Phigl=zeros(1,Maxorderl-1);
        
        dtheta=2*pi/1000;
        thetaint=linspace(0,2*pi,1001);
        drgl=interp1(thdrgl*pi/180,drgl,thetaint,'linear','extrap');
        
        for i=2:Maxorderl
            
            fint=drgl.*cos(i*pi/180*thetaint);
            ak=sum((fint(1:end-1,:)+fint(2:end,:))/2*dtheta)/2/pi;
            fint=drgl.*sin(i*pi/180*thetaint);
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

tilt_thu=0;
tilt_thl=0;

hrc=hui+thrt*aui+hli+thrb*ali;
hgc=hgi+(Db/2-Drg/2-gr_exp-arm)*(thgt+thgb);
gcl=hgc-hrc;

Pu=0; 
Pd=0; 
Pi=0;
Vp=0;   
Ap=-0;
Fi=-mr*Ap;

%=====================prepare input matrix======================%
inputs.pistonSM.beta_p=beta_p;
inputs.pistonSM.off=off;

inputs.bore.Db=Db;
inputs.bore.Ab=Ab;
inputs.bore.Phib=Phib;

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
inputs.ring.idr2g=idr2g;
inputs.ring.idg2r=idg2r;
inputs.ring.isi=isi;

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
inputs.FEM.fref=fref;
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
    PressInput=dlmread('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/pressure.txt'); %user has to provide this
    Tpo=PressInput(:,1);
    Pmo=PressInput(:,2);

    [~,~,Kapfs]=fscal_n4g(Tpo,Pmo,Rr,Er,Izr,gap,Nbe,1000);
end

%=====================solver======================%
Newton_tol=1e-6; %Tolerance for Newton-Raphson algorithm convergence (suggested value: 1e-6)
ItMax=100; %Maximum number of Newton-Raphson algorithm iterations (suggested value: 100)

Ug1=zeros(12*Nbnod,1);
Ug1(7:12:end-5)=alp;

F1=Finitial_fs_n4g(Ug1,Er,Gr,Iyy,Izz,Jt,alp,Rr,Telt,Nbe,1000,Kapfs);

inputs.ring.F1=F1;

[Kg_nd]=Kg_ring_curved_nd_n4g(Er,Gr,Iyy,Jt,Rr,Telt,Nbe,arm,ali,alo,isi,fref,href);

inputs.FEM.Kg_nd=Kg_nd;

UgIC_nd=zeros(12*Nbnod,1);
[ybn]=bore_dist_static_n4g(inputs);
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

UgIC_nd(1:12:12*Nbnod-11)=ybn(id1)+(Db/2-Rr-arm-sigmap)/uref;
UgIC_nd(2:12:12*Nbnod-10)=ybn(id2);
UgIC_nd(3:12:12*Nbnod-9)=ybn(id3);

UgIC_nd(4:12:12*Nbnod-8)=(4*omega*sigmagb-gcl/2)/href*ones(Nbnod,1);


if Pu>(Pi+Pd)
    beta=beta_p*cos(Tnod')+tilt_thu+1e-6; 
    UgIC_nd(7:12:12*Nbnod-5)=beta/href*Rr-sign(F1(7))*1e-6*sin((Tgg:Telt:2*pi-Tgg)/2)'/href*Rr;
    UgIC_nd(8:12:12*Nbnod-4)=-beta_p*sin(Tnod')/href*Rr*Telt-sign(F1(7))*1e-6*cos((Tgg:Telt:2*pi-Tgg)'/2)/href*Rr*Telt;
    UgIC_nd(9:12:12*Nbnod-3)=-beta_p*cos(Tnod')/href*Rr*Telt^2+sign(F1(7))*1e-6*sin((Tgg:Telt:2*pi-Tgg)'/2)/href*Rr*Telt^2;

else
    UgIC_nd(7:12:12*Nbnod-5)=-sign(F1(7))*1e-4*sin((Tgg:Telt:2*pi-Tgg)/2)/href*Rr;
    UgIC_nd(8:12:12*Nbnod-4)=-sign(F1(7))*1e-4*cos((Tgg:Telt:2*pi-Tgg)/2)/href*Rr*Telt;
    UgIC_nd(9:12:12*Nbnod-3)=sign(F1(7))*1e-4*sin((Tgg:Telt:2*pi-Tgg)/2)/href*Rr*Telt^2;
end

Ug_nd=UgIC_nd;

judge1=1;
judge2=1;
[NLS_0,yr,~,alr,hmin,fl,~,~,~,~]=Static_Error_n4g(Ug_nd,F1,inputs);
NLS_k=NLS_0;
j=0;
while j<=ItMax&&(judge1||judge2)
    j=j+1;
    [Ja_stac]=Jacobian_Estac_n4g(Ug_nd,inputs);
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
            [NLS_k,yr,~,alr,hmin,fl,~,~,~,~]=Static_Error_n4g(Ug_nd,F1,inputs);
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
            [NLS_k,yr,~,alr,hmin,fl,~,~,~,~]=Static_Error_n4g(Ug_nd,F1,inputs);
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

 N1=1-10*x.^3+15*x.^4-6*x.^5;
 N2=x-6*x.^3+8*x.^4-3*x.^5;
 N3=(x.^2)/2-3*(x.^3)/2+3*(x.^4)/2-(x.^5)/2;
 N4=10*x.^3-15*x.^4+6*x.^5;
 N5=-4*x.^3+7*x.^4-3*x.^5;
 N6=(x.^3)/2-x.^4+x.^5/2;
 N=[N1;N2;N3;N4;N5;N6];
 
 [ybn]=bore_dist_static_n4g(inputs);
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

 [zgnu,zgnl]=groove_dist_static_n4g(inputs);
 
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
 loc=zeros(1,Nr+1);
 lic=zeros(1,Nr+1);

 for j=1:Nbe
     ind=1+(j-1)*Npe:1+j*Npe;
     
     Uey=Ug_nd([12*(j-1)+1:12*(j-1)+3,12*j+1:12*j+3]);
     Uez=Ug_nd([12*(j-1)+4:12*(j-1)+6,12*j+4:12*j+6]);
     Uet=Ug_nd([12*(j-1)+7:12*(j-1)+9,12*j+7:12*j+9]);
     
     yrl=(Uey')*N*uref;
     zrl=(Uez')*N*href;
     alrl=(Uet')*N*href/Rr;
     beta=beta_p*cos((Tnod(j)+Telt*x));
     beta_thl=tilt_thl*ones(size(x)); 
     dz_g=(Rr*cos((Tnod(j)+Telt*x))-off)*beta_p;
     
     id=3*(j-1)+1:3*(j+1);
     yb=(ybn(id))'*N*uref;
     
     ybr(ind)=yb;
     yring(ind)=yrl;
     
     zgl=(zgnl(id))'*N*href; 
     
     Nby=50;

     yl=linspace(y1l,y2l,Nby+1);
     
     z_2l=repmat(dz_g+zgl-zrl,Nby+1,1);
     al_2l=repmat(beta+beta_thl-alrl,Nby+1,1);
 
     yyl=repmat(yl',1,Npe+1);
     
     hgl=gcl/2-z_2l-al_2l.*yyl;
     loc(ind)=hgl(end,:);
     lic(ind)=hgl(1,:);
 
 end
 
  display(['Ft(CBM) = ', num2str(sum(fl*Rr/Nr)),'N'])  %tengential load
 thee=linspace(Tgg*180/pi,360-Tgg*180/pi,Nr+1)';

 hmin=hmin';
 hmin=1e6*max(0,hmin);
 yr=1e6*yr';
 alr=180/pi*alr';
 loc=loc';
 loc=1e6*max(0,loc);
 lic=lic';
 lic=1e6*max(0,lic);
 
 figure(1)
 set(gcf,'Position',[314    58   703   916])
 subplot(311)
 hold on
 plot(thee(1:Npe:end),yr(1:Npe:end),'ro')
 plot(thee,yr,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Radial coordinate (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 subplot(312)
 hold on
 plot(thee(1:Npe:end),(hmin(1:Npe:end)),'ro')
 plot(thee,hmin,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Ring Liner Clearance (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 subplot(313)
 hold on
 plot(thee(1:Npe:end),alr(1:Npe:end),'ro')
 plot(thee,alr,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Ring Twist Angle (degree)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);

 figure(2)
 clf
 set(gcf,'Position',[144 64  1302 916])
 subplot(211)
 hold on
 plot(thee(1:Npe:end),loc(1:Npe:end),'ro')
 plot(thee,loc,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Lower OD Clearance (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 subplot(212)
 hold on
 plot(thee(1:Npe:end),lic(1:Npe:end),'ro')
 plot(thee,lic,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Lower ID Clearance (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);

 sigmagt=inputs.piston.sigmag1t;
 sigmagb=inputs.piston.sigmag1b;

fcyglob=zeros(1,Nr+1);
fczglob=zeros(1,Nr+1);
mcglob=zeros(1,Nr+1);
fclglob=zeros(1,Nr+1);
 
for j=1:Nbe
    
    ind=(Npe*(j-1)+1):(Npe*j+1);
    
    Uey=Ug_nd([12*(j-1)+1:12*(j-1)+3,12*j+1:12*j+3]);
    Uez=Ug_nd([12*(j-1)+4:12*(j-1)+6,12*j+4:12*j+6]);
    Uet=Ug_nd([12*(j-1)+7:12*(j-1)+9,12*j+7:12*j+9]);
    
    yre=(Uey')*N*uref;
    zre=(Uez')*N*href;
    
    alre=(Uet')*N*href/Rr;

    beta=beta_p*cos((Tnod(j)+Telt*x));

    beta_thu=tilt_thu*ones(size(x)); 

    beta_thl=tilt_thl*ones(size(x)); 
    dz_g=(Rr*cos((Tnod(j)+Telt*x))-off)*beta_p;
    
    id=3*(j-1)+1:3*(j+1);
    yb=(ybn(id))'*N*uref;
    
    zgu=(zgnu(id))'*N*href; 
    zgl=(zgnl(id))'*N*href;    
    
    [fcgu,mcgu]=ring_groove_contact_n4g(gcl,dz_g,zgu,zre,alre,alp,beta,(beta_thu+thgt-thrt),y1u,y2u,sigmagt,z,Pk_rg,omega,1,thrt,thrb);
    [fcgl,mcgl]=ring_groove_contact_n4g(gcl,dz_g,zgl,zre,alre,alp,beta,(beta_thl-thgb+thrb),y1l,y2l,sigmagb,z,Pk_rg,omega,-1,thrt,thrb);
    [fcl,mcl,~]=ring_liner_contact_n4g(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yre,alre,alp,PR,sigmap,z,Pk_rl,omega);
   
    Uef=Ug_nd([12*(j-1)+10:12*(j-1)+12,12*j+10:12*j+12]);

    ffc=(Uef')*N*fref;
    
    fcz=ffc+fcgu+fcgl;
    fczglob(ind)=fcz;

    fcy=-fcgu*thrt+fcgl*thrb+fcl;

    fclglob(ind)=fcl;
    fcyglob(ind)=fcy;
    mca=mcgu+fcgu*hui*thrt+mcgl+fcgl*hli*thrb+mcl+ffc*0;
    
    mcglob(ind)=mca;
    
end
 
 fy=fcyglob';
 fz=fczglob';
 m=mcglob';

figure(3)
set(gcf,'Position',[144 64  1302 916])
subplot(311)
hold on
plot(thee(1:Npe:end),-fy(1:Npe:end),'ro')
plot(thee,-fy,'linewidth',1.5)
xlabel('Circumferential Direction (degree)','FontSize',12)
ylabel('Radial Force (N/m)','FontSize',12)
set(gca,'XTick',[0:90:360])
xlim([0 360])
hold off
grid
set(gca,'FontSize',16);

subplot(312)
hold on
plot(thee(1:Npe:end),fz(1:Npe:end),'ro')
plot(thee,fz,'linewidth',1.5)
xlabel('Circumferential Direction (degree)','FontSize',12)
ylabel('Axial Force (N/m)','FontSize',12)
set(gca,'XTick',[0:90:360])
xlim([0 360])
hold off
grid
set(gca,'FontSize',16);

subplot(313)
hold on
plot(thee(1:Npe:end),m(1:Npe:end),'ro')
plot(thee,m,'linewidth',1.5)
xlabel('Circumferential Direction (degree)','FontSize',12)
ylabel('Twist moment (N)','FontSize',12)
set(gca,'XTick',[0:90:360])
xlim([0 360])
hold off
grid
set(gca,'FontSize',16);

mat_result=[thee,yr,hmin,alr,loc,lic,-fy,fz,m];
dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/stat_twist_result.txt',mat_result,'delimiter','\t','precision','%12.8f','newline','pc');  
end
