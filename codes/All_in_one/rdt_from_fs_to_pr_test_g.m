function rdt_from_fs_to_pr_test_g()
close all;

%this code computes the pressure distribution in circular bore
%the results are given within the variable res_force
%1st column contains angles in degrees and second one force per
%circumference length in N/m 

%the user either provides the free shape coordinates with first column
%containing angles in degrees and the second one containing the freeshape
%in millimeters
%or provides the free shape curvature with first column containing angles
%in degrees and the second one containing the freeshape curvature in 1/m

addpath(genpath(pwd))

%========================bore input==============================%
Db=95.25e-3;  %bore diameter, m
IsBDist=0; %1 to take into account bore thermal distortion, 0 otherwise
if IsBDist %considering bore distortion
    Ampdata=1; %1 to provide the amplitudes and phases directly or 0 to provide dr distribution
    if Ampdata
        Ab=[91.65 18 11.93 5.07]*1e-6;    %magnitude(m): 0th order, 2nd order, 3rd order, 4th order, etc
        Phib=[0.5642649 1.831025 1.1599458];  %phase(rad): 2nd order, 3rd order, 4th order, etc
    else
        
        prompt='Highest order in the Fourier series: ';
        Maxorder=input(prompt);
        
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
            
        dtheta=2*pi/10000;
        thetaint=linspace(0,2*pi,10001);
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
        
else %not considering bore distortion
    Ab=zeros(1,2); Phib=0;
end

%============================ring input=======================%
auo=2e-3;   %upper OD width, m
aui=2e-3;   %upper ID width, m
alo=2e-3;   %lower OD width, m
ali=2e-3;   %lower ID width, m
huo=1e-3;   %upper OD height, m
hui=1e-3;  %upper ID height, m
hlo=1e-3;  %lower OD height, m
hli=1e-3;  %lower ID height, m
thrt=0;           %ring upper flank angle, rad
thrb=0;           %ring lower flank angle, rad

rb1=0.2e-3;    %the lower edge width, m
rb2=0.7e-3;     %the upper edge width, m
rbn=-0;  %running face minimum point axial location, m
a10=0; 
a11=0; %Linear cofficient for lower edge shape factor:
a12=53;        %the lower edge shape factor, 1/m
a20=0;
a21=0; %Linear cofficient for upper edge shape factor:
a22=148;       %the upper edge shape factor, 1/m

gap=0.314e-3; %ring gap size when closed inside bore, m
arm=2e-3; %running face minimum point width, m

Ac=4e-6;   %area of cross-section, m^2
Izz=10.6667e-12; %principal moment of inertial in plane Izp, m^4
Iyy=2.6667e-12; %principal moment of inertial out of plane Iyp, m^4
Izr=10.6667e-12; %moment of inertial in plane Iz, m^4
Ip=Izz+Iyy;
Jt=7.328e-12;   %torsional factor, m^4
alp=0/180*pi;  %principal angle, rad

rhor=7820;       %ring density, kg/m^3
Er=100e9;        %ring Young's modulus, Pa
nur=0.3;         %ring Poisson Ratio
Gr=Er/(2*(1+nur));
alT=1.1e-5;       %thermal expansion coefficient, 1/K

Tg=0/180*pi;    %ring gap position, rad
rTemp=25;        %ring temperature, (C)

Rr=Db/2-arm;
Rr=Rr*(1+alT*(rTemp-25));

%=========================FEM input============================%
Nbe=16; %number of elements
Nbnod=Nbe+1;

Tgg=gap/2/Rr;

Telt=(2*pi-2*Tgg)/(Nbe);
Npe=1000; %number of points within a element for contact calculation
Nr=Nbe*Npe;

uref=5e-6;
href=20e-6;

Tnod=(Tgg:(Telt):2*pi-Tgg)+Tg;

Le=Rr*Telt;
Vnod=[Rr*cos(Tnod),Rr*sin(Tnod)];

mr=rhor*2*pi*Rr*Ac;

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
sigmagt=0.4e-6;            %ring groove upper flank roughness, m
sigmagb=0.4e-6;            %ring groove lower flank roughness, m
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
Pk_rl=cfct*K*A*2./((1-nul^2)/El+(1-nul^2)/Er);
Pk_rg=cfct*K*A*2./((1-nug^2)/Eg+(1-nur^2)/Er);

y1u=-aui;
y2u=arm-max((arm-auo),(Db-Drldu)/2+Ab(1)-exp_land_u);
lfu=y2u-y1u;
y1l=-ali;
y2l=arm-max((arm-alo),(Db-Drldl)/2+Ab(1)-exp_land_l);
lfl=y2l-y1l;
        
%======================groove input================%
IsGDist=0; % 1 to take into account groove thermal distortion, 0 otherwise:
if IsGDist %considering groove distortion
    Ampdatau=1; %1 to provide the amplitudes and phases directly or 0 to provide dz distribution
    if Ampdatau
        Agu=[20,0,0]*1e-6*0;   %2nd ring groove upper flank thermal distorsion magnitude 2nd, 3rd, 4th order, etc
        Phigu=[0,0,0]; %phase(rad): 2nd order, 3rd order, 4th order, etc
    else
        prompt='Number of orders starting from the 2nd one in the Fourier series for the upper groove deformation: ';
        Maxorderu=input(prompt);
        
        thdrgu=linspace(0,360,1e4)'; %ask for file containing first column of angles in degrees and second 
        dzgu=zeros(1,1e4)'; %column containing dz in microns
        
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
    Ampdatal=1; %1 to provide the amplitudes and phases directly or 0 to provide dz distribution
    if Ampdatal
        Agl=[20,0,0]*1e-6*0;   %2nd ring groove lower flank thermal distorsion magnitude 2nd, 3rd, 4th order, etc
        Phigl=[0,0,0]; %phase(rad): 2nd order, 3rd order, 4th order, etc
    else
        prompt='Number of orders starting from the 2nd one in the Fourier series for the lower groove deformation: ';
        Maxorderl=input(prompt);
        
        thdrgl=linspace(0,360,1e4)'; %ask for file containing first column of angles in degrees and second 
        dzgl=zeros(1,1e4)'; %column containing dz in microns
        
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

%=====================prepare input matrix======================%
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

%=============freeshape================%
IsFreeshape=0; %1 if user provides free shape coordinates, 0 if provides free shape curvature
if IsFreeshape %considering free shape coordinates

    %ask for free shape data with first column containing angles in degrees and
    %the second one containing the freeshape in millimeters

    fs=dlmread('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/fs_from_theo_pr_U.txt');

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
    
    %ask for free shape curvature with first column containing angles in degrees and
    %the second one containing the freeshape curvature in 1/m
    Kapfs=dlmread('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/kapfs_from_ov_from_fs_from_force_U_gap.txt');
    
    Tri=Kapfs(:,1)*pi/180;
    Kapfs=interp1(Tri,Kapfs(:,2),linspace(Tgg,2*pi-Tgg,Npe*Nbe+1)','linear','extrap');
end

%=====================solver======================%
Newton_tol=1e-6; %Tolerance for Newton-Raphson algorithm convergence (suggested value: 1e-6)
ItMax=100; %Maximum number of Newton-Raphson algorithm iterations (suggested value: 100)

Ug1=zeros(8*Nbnod,1);
Ug1(7:8:end-1)=alp;

F1=Finitial_fs_o(Ug1,Er,Gr,Iyy,Izz,Jt,alp,Rr,Telt,Nbe,Npe,Kapfs);

inputs.ring.F1=F1;

[Kg_nd]=Kg_ring_curved_nd_o(Er,Gr,Iyy,Jt,Rr,Telt,Nbe);

inputs.FEM.Kg_nd=Kg_nd;

UgIC_nd=zeros(8*Nbnod,1);
[ybn]=bore_dist_static_o(inputs);
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
[zgnu,~]=groove_dist_static_o(inputs);
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

judge1=1;
judge2=1;
NLS_0=Static_Error_o(Ug_nd,F1,inputs);
NLS_k=NLS_0;
j=0;

tic
while j<=ItMax&&(judge1||judge2)
    j=j+1;
    [Ja_stac]=Jacobian_Estac_o(Ug_nd,inputs);
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
            NLS_k=Static_Error_o(Ug_nd,F1,inputs);
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
            NLS_k=Static_Error_o(Ug_nd,F1,inputs);
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

Nt1=1-3*x.^2+2*x.^3;
Nt2=x-2*x.^2+x.^3;
Nt3=3*x.^2-2*x.^3;
Nt4=-x.^2+x.^3;

Nt=[Nt1;Nt2;Nt3;Nt4];

[ybn]=bore_dist_static_o(inputs);
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
res_force=[thee' fclr'];
dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/recovered_pressure_U_gap.txt',res_force,'delimiter','\t','precision','%12.8f','newline','pc');

figure(1)
clf
hold on
plot(thee,fclr,'Linewidth',1.5)
plot((Tnod-Tg)/pi*180,fclr(1:Npe:end),'ro')
xlabel('Circumferential Direction (degree)')
ylabel('Ring Liner Force (N/m)','FontSize',12)
set(gca,'XTick',[0:90:360])
hold off
grid on
set(gca,'FontSize',16);

end