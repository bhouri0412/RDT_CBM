function [rforce,rov]=rdt_from_fs_to_pr_and_ov(inputs)

%this code computes the radial linear force in ring fixed reference: 1st
%column contains angles in degrees and second one force per circumferential
%length in N/m. it also computes the ovality in usual representation ring 
%fixed reference: 1st column contains angles in degrees and second one 
%the coordintes in microns in usual representation 
%as inputs it needs a text file: 1st column contains angles in degrees in
%ring fixed reference and second one contains the free-shape coordinates
%in millimeters or its curvature in 1/m in usual representation

addpath(genpath(pwd))

%========================bore input==============================%
prompt='Bore diameter (mm): ';
Db=input(prompt)*1e-3;

prompt='Type 1 to take into account bore thermal distortion, 0 otherwise: ';
IsBDist=input(prompt);
if IsBDist
    
    prompt='Type 1 to provide the amplitudes and phases directly or 0 to provide dr distribution: ';
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
        
        prompt='Highest order in the Fourier seriesfor bore distortion: ';
        Norder=input(prompt);
        
        %ask for file containing bore distortion: first column contains
        %angles in degrees in bore fixed reference from thrust to thrust 
        %side and second column contains dr in microns
        dr=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\bore_dr.txt');  %user has to provide this
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
            
        dtheta=2*pi/10000;
        thetaint=linspace(0,2*pi,10001);
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

%============================ring input=======================%
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

prompt='Gap size when ring is closed within the bore (mm): ';
gap=input(prompt)*1e-3;

prompt='Gap size when ring is closed to ovality (under constant pressure) (mm): ';
gap2=input(prompt)*1e-3;

arm=inputs.ring.arm;
Ac=inputs.ring.Ac;
Izz=inputs.ring.Izz;
Iyy=inputs.ring.Iyy;
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

prompt='Ring gap position (deg): ';
Tg=input(prompt)*pi/180; 

prompt='Ring temperatrue (C): ';
rTemp=input(prompt);

Rr=Db/2-arm;
Rr=Rr*(1+alT*(rTemp-25));   

%=========================FEM input============================%
prompt='Number of elements for the FEM: ';
Nbe=input(prompt); 
Nbnod=Nbe+1; 
Tgg=gap/2/Rr;
Telt=(2*pi-2*Tgg)/(Nbe);
prompt='Number of points within one element: ';
Npe=input(prompt); 
Nr=Nbe*Npe;

uref=5e-6;  
href=20e-6;

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
prompt='Type 1 to take into account groove thermal distortion, 0 otherwise: ';
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

prompt='Type 1 to consider the linear force in computing the ovality equal to the average of the linear force in circular shape, 0 to provide its value: ';
IsPovAvg=input(prompt);
if ~IsPovAvg
    prompt='Constant linear force (N/m) to consider in computing the ovality (Pressure times the ring height (axial width))';
    Pov=input(prompt);
end

%=====================prepare local input matrix======================%
inputsl1.bore.Db=Db;
inputsl1.bore.Ab=Ab;
inputsl1.bore.Phib=Phib;
inputsl1.ring.auo=auo;
inputsl1.ring.aui=aui;
inputsl1.ring.alo=alo;
inputsl1.ring.ali=ali;
inputsl1.ring.huo=huo;
inputsl1.ring.hui=hui;
inputsl1.ring.hlo=hlo;
inputsl1.ring.hli=hli;
inputsl1.ring.thrt=thrt;           
inputsl1.ring.thrb=thrb;           
inputsl1.ring.rb1=rb1;
inputsl1.ring.rb2=rb2;
inputsl1.ring.rbn=rbn;
inputsl1.ring.a10=a10;
inputsl1.ring.a11=a11;
inputsl1.ring.a12=a12;
inputsl1.ring.a20=a20;
inputsl1.ring.a21=a21;
inputsl1.ring.a22=a22;
inputsl1.ring.arm=arm;
inputsl1.ring.gap=gap;
inputsl1.ring.rhor=rhor;
inputsl1.ring.Er=Er;
inputsl1.ring.alT=alT;
inputsl1.ring.nur=nur;
inputsl1.ring.Gr=Gr;
inputsl1.ring.Ac=Ac;
inputsl1.ring.Izz=Izz;
inputsl1.ring.Iyy=Iyy;
inputsl1.ring.Izr=Izr;
inputsl1.ring.Ip=Ip;
inputsl1.ring.Jt=Jt;
inputsl1.ring.alp=alp;
inputsl1.ring.Rr=Rr;
inputsl1.ring.Le=Le;
inputsl1.ring.mr=mr;
inputsl1.ring.Tg=Tg;
inputsl1.ring.rTemp=rTemp;

inputsl1.piston.Drldu=Drldu;
inputsl1.piston.Drg=Drg;
inputsl1.piston.Drldl=Drldl;
inputsl1.piston.hgi=hgi;
inputsl1.piston.thgt=thgt;
inputsl1.piston.thgb=thgb;
inputsl1.piston.sigmag1t=sigmagt;
inputsl1.piston.sigmag1b=sigmagb;
inputsl1.piston.nug=nug;
inputsl1.piston.Eg=Eg;
inputsl1.piston.exp_land_u=exp_land_u;
inputsl1.piston.exp_land_l=exp_land_l;
inputsl1.piston.gr_exp=gr_exp;

inputsl1.liner.El=El;
inputsl1.liner.nul=nul;
inputsl1.liner.sigmap=sigmap;
inputsl1.liner.PR=PR;

inputsl1.FEM.Nbe=Nbe;
inputsl1.FEM.Nbnod=Nbnod;
inputsl1.FEM.Telt=Telt;
inputsl1.FEM.Npe=Npe;
inputsl1.FEM.Nr=Nr;
inputsl1.FEM.Le=Le;
inputsl1.FEM.Tnod=Tnod;
inputsl1.FEM.Vnod=Vnod;
inputsl1.FEM.uref=uref;
inputsl1.FEM.href=href;

inputsl1.contact.z=z;
inputsl1.contact.Pk_rl=Pk_rl;
inputsl1.contact.Pk_rg=Pk_rg;
inputsl1.contact.omega=omega;
inputsl1.contact.fc_dry=fc_dry;
inputsl1.contact.cfct=cfct;

inputsl1.groove.y1u=y1u;
inputsl1.groove.y2u=y2u;
inputsl1.groove.lfu=lfu;
inputsl1.groove.y1l=y1l;
inputsl1.groove.y2l=y2l;
inputsl1.groove.lfl=lfl;
inputsl1.groove.gcl=gcl;
inputsl1.groove.Agu=Agu;
inputsl1.groove.Phigu=Phigu;
inputsl1.groove.Agl=Agl;
inputsl1.groove.Phigl=Phigl;
inputsl1.groove.tilt_thu=tilt_thu;
inputsl1.groove.tilt_thl=tilt_thl;

%=============freeshape================%
prompt='Type 1 to provide the free shape coordinates and 0 to provide the free shape curvature: ';
IsFreeshape=input(prompt);
if IsFreeshape

    %ask for file containing free-shape coordinates in ring fixed
    %reference: first column containins angles in degrees and second 
    %one contains the freeshape coordinates in usual representation in
    %millimeters
    fs=dlmread('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\Texts\fs_from_theo_pr_TC.txt'); %user has to provide this

    Tri=fs(:,1)*pi/180;
    rfsim=fs(:,2)*1e-3;
    Inod=1:Npe:Npe*Nbe+1;
    dTrfs=(Tri(end)-Tri(1))/(Npe*Nbe);
    Trfs=(Tri(1):dTrfs:Tri(length(Tri)))';
    Teltfs=(Trfs(end)-Trfs(1))/Nbe;
    Tnodfs=Trfs(Inod);
    rfsi=interp1(Tri,rfsim,Trfs);
    [Ccs]=Cs_3rd_cont_o(Tnodfs);
    [rnodfs]=ls_rfsi_o(rfsi,Nbnod,Teltfs,Ccs);
    [Kapfs,~]=Kfs_cal_o(Trfs,rfsi,rnodfs,Nbe,Npe,Teltfs);
else
    
    %ask for file containing free-shape curvature in ring fixed
    %reference: first column containins angles in degrees and second 
    %one contains the freeshape curvature usual representation in 1/m
    Kapfs=dlmread('C:\Users\YangLiu\Dropbox (MIT)\Aziz Perso\papers mit\RA\Thesis\Texts\Kapfs_from_theo_pr_TC.txt');  %user has to provide this
    Tri=Kapfs(:,1)*pi/180;
    Kapfs=interp1(Tri,Kapfs(:,2),linspace(Tgg,2*pi-Tgg,Npe*Nbe+1)','linear','extrap');
end

%=====================solver to compute the pressure======================%
Ug1=zeros(8*Nbnod,1);
Ug1(7:8:end-1)=alp;

F1=Finitial_fs_o(Ug1,Er,Gr,Iyy,Izz,Jt,alp,Rr,Telt,Nbe,Npe,Kapfs);

inputsl1.ring.F1=F1;

[Kg_nd]=Kg_ring_curved_nd_o(Er,Gr,Iyy,Jt,Rr,Telt,Nbe);

inputsl1.FEM.Kg_nd=Kg_nd;

UgIC_nd=zeros(8*Nbnod,1);
[ybn]=bore_dist_static_o(inputsl1);
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
[zgnu,~]=groove_dist_static_o(inputsl1);
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

prompt='Tolerance for Newton-Raphson algorithm convergence used to compute the ovality (suggested value: 1e-6): ';
Newton_tol2=input(prompt); 
prompt='Maximum number of Newton-Raphson algorithm iterations used to compute the ovality (suggested value: 100): ';
ItMax2=input(prompt);

prompt='Tolerance for Newton-Raphson algorithm convergence used to compute the force distribution (suggested value: 1e-6): ';
Newton_tol=input(prompt);
prompt='Maximum number of Newton-Raphson algorithm iterations used to compute the force distribution (suggested value: 100): ';
ItMax=input(prompt); 
judge1=1;
judge2=1;
NLS_0=Static_Error_o(Ug_nd,F1,inputsl1);
NLS_k=NLS_0;
j=0;

while j<=ItMax&&(judge1||judge2)
    j=j+1;
    [Ja_stac]=Jacobian_Estac_o(Ug_nd,inputsl1);
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
            NLS_k=Static_Error_o(Ug_nd,F1,inputsl1);
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
            NLS_k=Static_Error_o(Ug_nd,F1,inputsl1);
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
hold on
plot(thee,fclr,'Linewidth',1.5)
plot((Tnod-Tg)/pi*180,fclr(1:Npe:end),'ro')
xlabel('Circumferential Direction (degree)')
ylabel('Ring Liner Force (N/m)','FontSize',12)
set(gca,'XTick',[0:90:360])
hold off
grid on
set(gca,'FontSize',16);

if IsPovAvg
    Pov=mean(fclr);
end

Tgg2=gap2/2/Rr;
Telt2=(2*pi-2*Tgg2)/(Nbe);
Kapfs=interp1((Tgg:Telt:2*pi-Tgg)',Kapfs,(Tgg2:Telt2:2*pi-Tgg2)','linear','extrap');
Tgg=Tgg2;
Telt=Telt2;
Le=Rr*Telt;

Ug=zeros(3*(Nbe+1),1);

%=====================solver to compute the ovality======================% 
judge1=1;
NLS_0=NewtonOV_Error(Ug,Er,Izr,Rr,Nbe,Npe,Le,Kapfs',Pov);
NLS_k=NLS_0;
j=0;

while j<=ItMax2&&(judge1)
    j=j+1;
    [Ja_ov]=Jacobian_ov(Nbe,Npe,Le,Er,Izr,Rr,Ug,Kapfs');
   
    Nstep=-Ja_ov\NLS_k; 
    Ug=Ug+Nstep;

    NLS_k=NewtonOV_Error(Ug,Er,Izr,Rr,Nbe,Npe,Le,Kapfs',Pov);
    judge1=norm(NLS_k)>Newton_tol2;

    norm_NLS_k_ov(j)=norm(NLS_k);
    
    figure(1)
    cla
    semilogy(norm_NLS_k_ov)
    pause(0.0001)
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
hold on
plot(Tr',ov'*1e6,'LineWidth',1.5)
xlabel('Circumferential Direction (degree)')
ylabel('Ovality (\mu m)')
title('Ovality in usual representation');
set(gca,'XTick',[0:90:360])
set(gca,'FontSize',16);
grid on

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