function [res_thee,res_yr,res_hmin,res_alr,res_loc,res_lic,res_fy,res_fz,res_m]=rdt_static_twist(inputs)

%this code performs the static twist calculations under FixOD or FixID Constraint

% it outputs the node's circumferential location (in microns), ring-liner
% clearance (microns), twist angle (degree); lower OD clearance (microns),
% lower ID clearance (microns)
%radial force (N/m), axial force (N/m) and twist moment (N) 

addpath(genpath(pwd))


IsBDist=0;

IsPTilt=0;

IsGDef=0;

%=============features================================%
prompt='Type 1 to provide free shape or its curvature and 0 to provide radial pressure distribution for the circular shape: ';
IsFreeshape=input(prompt);
inputsl5.features.IsBDist=IsBDist;
inputsl5.features.IsPTilt=IsPTilt;
inputsl5.features.IsGDef=IsGDef;
inputsl5.features.IsFreeshape=IsFreeshape;

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
Izr=inputs.ring.Izr; 
Ip=inputs.ring.Ip; 
Jt=inputs.ring.Jt; 
alp=inputs.ring.alp; 
isi=0*(alp<=0)+1*(alp>0);

prompt='Ring density (kg/m^3): ';
rhor=input(prompt); 

prompt='Ring Young modulus (GPa): ';
Er=input(prompt)*1e9;

prompt='Ring Poisson Ratio: ';
nur=input(prompt);

Gr=Er/(2*(1+nur));

alT=1.1e-5; 

Tg=0;

rTemp=25;

Rro=Db/2-arm;
Rr=Rro*(1+alT*(rTemp-25));    

%=========================FEM input============================%
prompt='Number of elements for the FEM: ';
Nbe=input(prompt); 
Nbnod=Nbe+1; 
Tgg=gap/2/Rr;
Telt=(2*pi-2*Tgg)/(Nbe);
prompt='Number of points within one element for contact and lubrication force calculation: ';
Npe=input(prompt); 
Nr=Nbe*Npe; 
Nge=Npe; 
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
    
    prompt='Type 1 to provide the amplitudes and phases directly or 0 to provide dr distribution: ';
    Ampdata=input(prompt);
    if Ampdata
        prompt='Highest order in the Fourier series: ';
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
        
        prompt='Highest order in the Fourier series: ';
        Maxorder=input(prompt);
        
        thdr=linspace(0,360,1e4)'; 
        dr=zeros(1,1e4)'; 
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
    prompt='Piston tilt angle (deg): ';
    beta_p=input(prompt)*pi/180; 
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
prompt='Liner Young modulus (GPa): ';
El=input(prompt)*1e9; 

prompt='Liner Poisson Ratio: ';
nul=input(prompt);

prompt='Plateau Ratio: ';
PR=input(prompt);

prompt='Liner surface roughness standard deviation (\mu m): ';
sigmap=input(prompt)*1e-6; 

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
prompt='z coefficient for the asperity ring/liner and ring/lower plate contact interaction (6.804 is the adopted value in the simplified formulation): ';
z=input(prompt); 
prompt='K coefficient for the asperity ring/liner and ring/lower plate contact interaction (1.198e-4 is the adopted value in the simplified formulation): ';
K=input(prompt); 
prompt='A coefficient for the asperity ring/liner and ring/lower plate contact interaction (4.4068e-5 is the adopted value in the simplified formulation): ';
A=input(prompt);
prompt='Omega coefficient for the asperity ring/liner and ring/lower plate contact interaction (4 is the adopted value in the simplified formulation): ';
omega=input(prompt);
fc_dry=0.13;
cfct=1;
Pk_rl=cfct*K*A*2./((1-nul^2)/El+(1-nul^2)/Er);
Pk_rg=cfct*K*A*2./((1-nug^2)/Eg+(1-nur^2)/Er);

y1u=-aui;
y2u=arm-max((arm-auo),(Db-Drldu)/2+Ab(1)-exp_land_u);
lfu=y2u-y1u;
y1l=-ali;
y2l=arm-max((arm-alo),(Db-Drldl)/2+Ab(1)-exp_land_l);
lfl=y2l-y1l;

if IsGDef
    prompt='Type 1 to provide the amplitudes and phases directly for the upper groove deformation or 0 to provide the axial coordinate distribution: ';
    Ampdatau=input(prompt);
    if Ampdatau
        prompt='Number of orders starting from the 2nd one in the Fourier series for the upper groove deformation: ';
        Norder=input(prompt);
        Agu=zeros(1,Norder);
        Phigu=zeros(1,Norder);
        for i=1:Norder
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
        Maxorderu=input(prompt);
        
        thdrgu=linspace(0,360,1e4)';  
        dzgu=zeros(1,1e4)';
        
        dzgu=dzgu'*1e-6;
        thdrgu=thdrgu';
        
        Agu=zeros(1,Maxorderu);
        Phigu=zeros(1,Maxorderu);
        
        dtheta=2*pi/10000;
        thetaint=linspace(0,2*pi,10001);
        dzgu=interp1(thdrgu*pi/180,dzgu,thetaint,'linear','extrap');
        
        for i=1:Maxorderu
            
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
        Norder=input(prompt);
        Agl=zeros(1,Norder);
        Phigl=zeros(1,Norder);
        for i=1:Norder
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
        Maxorderl=input(prompt);
        
        thdrgl=linspace(0,360,1e4)'; 
        dzgl=zeros(1,1e4)'; 
        
        dzgl=dzgl'*1e-6;
        thdrgl=thdrgl';
        
        Agl=zeros(1,Maxorderl);
        Phigl=zeros(1,Maxorderl);
        
        dtheta=2*pi/10000;
        thetaint=linspace(0,2*pi,10001);
        dzgl=interp1(thdrgl*pi/180,dzgl,thetaint,'linear','extrap');
        
        for i=1:Maxorderl
            
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

%=====================prepare local input matrix======================%
inputsl5.pistonSM.beta_p=beta_p;
inputsl5.pistonSM.off=off;

inputsl5.bore.Db=Db;
inputsl5.bore.Ab=Ab;
inputsl5.bore.Phib=Phib;

inputsl5.ring.auo=auo;
inputsl5.ring.aui=aui;
inputsl5.ring.alo=alo;
inputsl5.ring.ali=ali;
inputsl5.ring.huo=huo;
inputsl5.ring.hui=hui;
inputsl5.ring.hlo=hlo;
inputsl5.ring.hli=hli;
inputsl5.ring.thrt=thrt;           
inputsl5.ring.thrb=thrb;           
inputsl5.ring.rb1=rb1;
inputsl5.ring.rb2=rb2;
inputsl5.ring.rbn=rbn;
inputsl5.ring.a10=a10;
inputsl5.ring.a11=a11;
inputsl5.ring.a12=a12;
inputsl5.ring.a20=a20;
inputsl5.ring.a21=a21;
inputsl5.ring.a22=a22;
inputsl5.ring.arm=arm;
inputsl5.ring.gap=gap;
inputsl5.ring.rhor=rhor;
inputsl5.ring.Er=Er;
inputsl5.ring.alT=alT;
inputsl5.ring.nur=nur;
inputsl5.ring.Gr=Gr;
inputsl5.ring.Ac=Ac;
inputsl5.ring.Izz=Izz;
inputsl5.ring.Iyy=Iyy;
inputsl5.ring.Izr=Izr;
inputsl5.ring.Ip=Ip;
inputsl5.ring.Jt=Jt;
inputsl5.ring.alp=alp;
inputsl5.ring.Rr=Rr;
inputsl5.ring.Le=Le;
inputsl5.ring.mr=mr;
inputsl5.ring.Tg=Tg;
inputsl5.ring.rTemp=rTemp;
inputsl5.ring.idr2g=idr2g;
inputsl5.ring.idg2r=idg2r;
inputsl5.ring.isi=isi;

inputsl5.piston.Drldu=Drldu;
inputsl5.piston.Drg=Drg;
inputsl5.piston.Drldl=Drldl;
inputsl5.piston.hgi=hgi;
inputsl5.piston.thgt=thgt;
inputsl5.piston.thgb=thgb;
inputsl5.piston.sigmag1t=sigmagt;
inputsl5.piston.sigmag1b=sigmagb;
inputsl5.piston.nug=nug;
inputsl5.piston.Eg=Eg;
inputsl5.piston.exp_land_u=exp_land_u;
inputsl5.piston.exp_land_l=exp_land_l;
inputsl5.piston.gr_exp=gr_exp;

inputsl5.liner.El=El;
inputsl5.liner.nul=nul;
inputsl5.liner.sigmap=sigmap;
inputsl5.liner.PR=PR;

inputsl5.FEM.Nbe=Nbe;
inputsl5.FEM.Nbnod=Nbnod;
inputsl5.FEM.Telt=Telt;
inputsl5.FEM.Npe=Npe;
inputsl5.FEM.Nr=Nr;
inputsl5.FEM.Le=Le;
inputsl5.FEM.Tnod=Tnod;
inputsl5.FEM.Vnod=Vnod;
inputsl5.FEM.uref=uref;
inputsl5.FEM.href=href;
inputsl5.FEM.fref=fref;
inputsl5.FEM.Nge=Nge;
inputsl5.FEM.Ng=Ng;

inputsl5.contact.z=z;
inputsl5.contact.Pk_rl=Pk_rl;
inputsl5.contact.Pk_rg=Pk_rg;
inputsl5.contact.omega=omega;
inputsl5.contact.fc_dry=fc_dry;
inputsl5.contact.cfct=cfct;

inputsl5.groove.y1u=y1u;
inputsl5.groove.y2u=y2u;
inputsl5.groove.lfu=lfu;
inputsl5.groove.y1l=y1l;
inputsl5.groove.y2l=y2l;
inputsl5.groove.lfl=lfl;
inputsl5.groove.gcl=gcl;
inputsl5.groove.Agu=Agu;
inputsl5.groove.Phigu=Phigu;
inputsl5.groove.Agl=Agl;
inputsl5.groove.Phigl=Phigl;
inputsl5.groove.tilt_thu=tilt_thu;
inputsl5.groove.tilt_thl=tilt_thl;

inputsl5.engine.Vp=Vp;
inputsl5.engine.Ap=Ap;
inputsl5.ring.Fi=Fi;
inputsl5.engine.Pu=Pu;
inputsl5.engine.Pd=Pd;
inputsl5.engine.Pi=Pi;

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
        [Ccs]=Cs_3rd_cont_n4g(Tnodfs);
        [rnodfs]=ls_rfsi_n4g(rfsi,Nbnod,Teltfs,Ccs);
        [Kapfs,~]=Kfs_cal_n4g(Trfs,rfsi,rnodfs,Nbe,1000,Teltfs);
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

    [~,~,Kapfs]=fscal_n4g(Tpo,Pmo,Rr,Er,Izr,gap,Nbe,1000);
end

%=====================solver======================%
Ug1=zeros(8*Nbnod,1);
Ug1(7:8:end-1)=alp;

F1=Finitial_fs_n4g(Ug1,Er,Gr,Iyy,Izz,Jt,alp,Rr,Telt,Nbe,1000,Kapfs);

inputsl5.ring.F1=F1;

[Kg_nd]=Kg_ring_curved_nd_n4g(Er,Gr,Iyy,Jt,Rr,Telt,Nbe,arm,ali,alo,isi,fref,href);

inputsl5.FEM.Kg_nd=Kg_nd;

UgIC_nd=zeros(12*Nbnod,1);
[ybn]=bore_dist_static_n4g(inputsl5);
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

prompt='Tolerance for Newton-Raphson algorithm convergence (suggested value: 1e-6): ';
Newton_tol=input(prompt); 
prompt='Maximum number of Newton-Raphson algorithm iterations (suggested value: 100): ';
ItMax=input(prompt);
judge1=1;
judge2=1;
[NLS_0,yr,~,alr,hmin,fl,~,~,~,~]=Static_Error_n4g(Ug_nd,F1,inputsl5);
NLS_k=NLS_0;
j=0;
while j<=ItMax&&(judge1||judge2)
    j=j+1;
    [Ja_stac]=Jacobian_Estac_n4g(Ug_nd,inputsl5);
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
            [NLS_k,yr,~,alr,hmin,fl,~,~,~,~]=Static_Error_n4g(Ug_nd,F1,inputsl5);
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
            [NLS_k,yr,~,alr,hmin,fl,~,~,~,~]=Static_Error_n4g(Ug_nd,F1,inputsl5);
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

 [ybn]=bore_dist_static_n4g(inputsl5);
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

 [zgnu,zgnl]=groove_dist_static_n4g(inputsl5);
 
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
 
  display(['Ft(CBM) = ', num2str(sum(fl*Rr/Nr)),'N']) %tengential load
 thee=linspace(Tgg*180/pi,360-Tgg*180/pi,Nr+1)';
 
 hmin=hmin';
 hmin=max(0,hmin);
 yr=yr';
 alr=alr';
 loc=loc';
 loc=max(0,loc);
 lic=lic';
 lic=max(0,lic);
 res_thee=thee;
 res_hmin=hmin*1e6;
 res_yr=yr*1e6;
 res_alr=alr/pi*180;
 res_loc=loc*1e6;
 res_lic=lic*1e6;
 
 figure(1)
 set(gcf,'Position',[314    58   703   916])
 subplot(311)
 hold on
 plot(thee(1:Npe:end),yr(1:Npe:end)*1e6,'ro')
 plot(thee,yr*1e6,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Radial coordinate (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 subplot(312)
 hold on
 plot(thee(1:Npe:end),(hmin(1:Npe:end))*1e6,'ro')
 plot(thee,hmin*1e6,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Ring Liner Clearance (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 subplot(313)
 hold on
 plot(thee(1:Npe:end),alr(1:Npe:end)/pi*180,'ro')
 plot(thee,alr/pi*180,'linewidth',1.5)
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
 plot(thee(1:Npe:end),loc(1:Npe:end)*1e6,'ro')
 plot(thee,loc*1e6,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Lower OD Clearance (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 subplot(212)
 hold on
 plot(thee(1:Npe:end),lic(1:Npe:end)*1e6,'ro')
 plot(thee,lic*1e6,'linewidth',1.5)
 xlabel('Circumferential Direction (degree)','FontSize',12)
 ylabel('Lower ID Clearance (\mum)','FontSize',12)
 set(gca,'XTick',[0:90:360])
 xlim([0 360])
 hold off
 grid
 set(gca,'FontSize',16);
 
 sigmagt=inputsl5.piston.sigmag1t;
 sigmagb=inputsl5.piston.sigmag1b;

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
 
 res_fy=-fy;
 res_fz=fz;
 res_m=m;

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
 
end