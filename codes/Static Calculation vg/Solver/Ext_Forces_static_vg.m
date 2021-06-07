function [Fext_nd,yr,zr,alr,hmin,fl,z0]=Ext_Forces_static_vg(Ug_nd,inputs)
%Calculation of nodal resultants of external forces acting on the ring (static
%analysis)

%-------------- Retrieving parameters from input structure ---------------%
Db=inputs.bore.Db;
temp=inputs.bore.templ;
dt=1;

%engine
Vp=inputs.engine.Vp;
Pu=inputs.engine.Pu;
Pd=inputs.engine.Pd;
Pi=inputs.engine.Pi;

%oil data
% zk=inputs.oil.zk;
% temp1=inputs.oil.temp1;
% temp2=inputs.oil.temp2;
% rho_oil=inputs.oil.rho_oil;
mu_oilv=inputs.oil.mu_oil;
%Meshing
Telt=inputs.FEM.Telt; %element size (rad)
Tnod=inputs.FEM.Tnod; %node angular position (rad)
Nbe=inputs.FEM.Nbe; %number of elements
Nbnod=inputs.FEM.Nbnod; %number of nodes
Nr=inputs.FEM.Nr;
Ng=inputs.FEM.Ng;
Npe=inputs.FEM.Npe; %number of contact points within an element
uref=inputs.FEM.uref;
href=inputs.FEM.href;
Nge=inputs.FEM.Nge;

%ring
Er=inputs.ring.Er;
alT=inputs.ring.alT;
Gr=inputs.ring.Gr;
Izz=inputs.ring.Izz;
Iyy=inputs.ring.Iyy;
Jt=inputs.ring.Jt;
Rr=inputs.ring.Rr; %ring neutral fiber radius
Le=inputs.ring.Le;
huo=inputs.ring.huo;
hlo=inputs.ring.hlo;
arm=inputs.ring.arm;
ali=inputs.ring.ali;
aui=inputs.ring.aui;
ar=arm+max(ali,aui);
rb1=inputs.ring.rb1;
rb2=inputs.ring.rb2;
rbn=inputs.ring.rbn;
a10=inputs.ring.a10;
a11=inputs.ring.a11;
a12=inputs.ring.a12;
a20=inputs.ring.a20;
a21=inputs.ring.a21;
a22=inputs.ring.a22;
dTemp=inputs.ring.dTemp;
hui=inputs.ring.hui;
hli=inputs.ring.hli;
thrt=inputs.ring.thrt;
thrb=inputs.ring.thrb;
Fi=inputs.ring.Fi;
Tg=inputs.ring.Tg;

%groove
gcl=inputs.groove.gcl;
y1u=inputs.groove.y1u;
y2u=inputs.groove.y2u;
y1l=inputs.groove.y1l;
y2l=inputs.groove.y2l;
lfu=inputs.groove.lfu;
lfl=inputs.groove.lfl;
tilt_thu=inputs.groove.tilt_thu; %groove thermal tilting - upper flank
tilt_thl=inputs.groove.tilt_thl; %groove thermal tilting - lower flank

%piston
sigmagt=inputs.piston.sigmag1t;
sigmagb=inputs.piston.sigmag1b;
%Drg=inputs.piston.Drg;
%hgi=inputs.piston.hgi;
%gr_exp=inputs.piston.gr_exp;
thgt=inputs.piston.thgt;
thgb=inputs.piston.thgb;
hr=hui+thrt*aui+hli+thrb*ali;     %at centroid
%hg=hgi+(Db/2-Drg/2-gr_exp-arm)*(thgt+thgb);   %at centroid
%hot=0e-6;
%hob=0e-6;

%liner
sigmap=inputs.liner.sigmap;
PR=inputs.liner.PR;
Ph_OCR=inputs.liner.Ph;
KOCR=inputs.liner.KOCR;
% muUo=inputs.liner.muUo;
% c1=inputs.liner.c1;
% c2=inputs.liner.c2;
% c3=inputs.liner.c3;
ap=inputs.liner.ap;
cp1=inputs.liner.cp1;
cp2=inputs.liner.cp2;
F0=inputs.liner.F0;
cf1=inputs.liner.cf1;
cf2=inputs.liner.cf2;
rw=inputs.liner.rw;
muUot=inputs.liner.muUot;
b=inputs.liner.b;
c=inputs.liner.c;
d=inputs.liner.d;
e=inputs.liner.e;

%contact
z=inputs.contact.z;
omega=inputs.contact.omega;
Pk_rg=inputs.contact.Pk_rg;
Pk_rl=inputs.contact.Pk_rl;
fc=inputs.contact.fc_dry;

%lubrication
hrls=inputs.lubrication.hrls;
isfuela=inputs.lubrication.isfuel;
hot=inputs.lubrication.hot;
hob=inputs.lubrication.hob;
minoil=inputs.lubrication.minoil;
oilthreshold=inputs.lubrication.oilthreshold;
viscosityfactor=inputs.lubrication.viscosityfactor;
MT=1/12*alT*Er*hr*ar^2*mean(dTemp);

%secondary motion
beta_p=inputs.pistonSM.beta_p; %piston tilt - current step
off=inputs.pistonSM.off; %pin offset
%-------------------------------------------------------------------------%
shiftstep=round(Tg/2/pi*Nr);
%temp_liner=circshift(temp,[0 -shiftstep]);
temp_liner=temp;

%shift to ring reference
%isfuel=circshift(isfuela,[0 -shiftstep]);
isfuel=isfuela;

hrls=hrls.*(hrls>=(minoil))+minoil*(hrls<(minoil));
hrl=2*hrls;
%hrl=circshift(2*hrls,[0 -shiftstep]);
Kp=cp1+cp2*hrl/sigmap;
Kf=cf1+cf2*hrl/sigmap;

shiftstepg=round(Tg/2/pi*Ng);
%Pu=circshift(Pu,[0 -shiftstepg]);
%Pd=circshift(Pd,[0 -shiftstepg]);
%Pi=circshift(Pi,[0 -shiftstepg]);
%hot=circshift(hot,[0 -shiftstepg]);
%hob=circshift(hob,[0 -shiftstepg]);
%mu_oilv=circshift(mu_oilv,[0 -shiftstepg]);

%discretization inside the element
dx=1/Npe;
x=0:dx:1;

dxg=1/Nge;
xg=0:dxg:1;

%Shape functions for z displacements
N1=1-10*x.^3+15*x.^4-6*x.^5;
N2=x-6*x.^3+8*x.^4-3*x.^5;
N3=(x.^2)/2-3*(x.^3)/2+3*(x.^4)/2-(x.^5)/2;
N4=10*x.^3-15*x.^4+6*x.^5;
N5=-4*x.^3+7*x.^4-3*x.^5;
N6=(x.^3)/2-x.^4+x.^5/2;
N=[N1;N2;N3;N4;N5;N6];
N_int=(N(:,1:end-1)+N(:,2:end))/2;

dN1=-30*x.^2+60*x.^3-30*x.^4;
dN2=1-18*x.^2+32*x.^3-15*x.^4;
dN3=x-9*(x.^2)/2+6*x.^3-5*(x.^4)/2;
dN4=30*x.^2-60*x.^3+30*x.^4;
dN5=-12*x.^2+28*x.^3-15*x.^4;
dN6=3*(x.^2)/2-4*x.^3+5*x.^4/2;
dN=[dN1;dN2;dN3;dN4;dN5;dN6];
dN_int=(dN(:,1:end-1)+dN(:,2:end))/2;

N1g=1-10*xg.^3+15*xg.^4-6*xg.^5;
N2g=xg-6*xg.^3+8*xg.^4-3*xg.^5;
N3g=(xg.^2)/2-3*(xg.^3)/2+3*(xg.^4)/2-(xg.^5)/2;
N4g=10*xg.^3-15*xg.^4+6*xg.^5;
N5g=-4*xg.^3+7*xg.^4-3*xg.^5;
N6g=(xg.^3)/2-xg.^4+xg.^5/2;
N_g=[N1g;N2g;N3g;N4g;N5g;N6g];
Ng_int=(N_g(:,1:end-1)+N_g(:,2:end))/2;

%Shape functions for twist displacements
Nt1=1-3*x.^2+2*x.^3;
Nt2=x-2*x.^2+x.^3;
Nt3=3*x.^2-2*x.^3;
Nt4=-x.^2+x.^3;
Nt=[Nt1;Nt2;Nt3;Nt4];
Nt_int=(Nt(:,1:end-1)+Nt(:,2:end))/2;

Nt1g=1-3*xg.^2+2*xg.^3;
Nt2g=xg-2*xg.^2+xg.^3;
Nt3g=3*xg.^2-2*xg.^3;
Nt4g=-xg.^2+xg.^3;
Ntg=[Nt1g;Nt2g;Nt3g;Nt4g];
Ntg_int=(Ntg(:,1:end-1)+Ntg(:,2:end))/2;

%Initializing vector of external forces
Fext_nd=zeros(8*Nbnod,1);

%Bore distortion at node loactions
[ybn]=bore_dist_static_vg(inputs);

Tnodd=(0:2*pi/Nbe:2*pi)+Tg;

%indexes for filling nodal displacement matrices
id1=1:3:3*Nbnod-2;
id2=2:3:3*Nbnod-1;
id3=3:3:3*Nbnod;

ybn1=ybn(id1);
ybn2=ybn(id2);
ybn3=ybn(id3);

ybn(id1)=interp1(Tnodd,ybn1,Tnod);
ybn(id2)=interp1(Tnodd,ybn2,Tnod);
ybn(id3)=interp1(Tnodd,ybn3,Tnod);

%Groove distortion at node locations
[zgnu,zgnl]=groove_dist_static_vg(inputs);

zgnu1=zgnu(id1);
zgnu2=zgnu(id2);
zgnu3=zgnu(id3);

zgnu(id1)=interp1(Tnodd,zgnu1,Tnod);
zgnu(id2)=interp1(Tnodd,zgnu2,Tnod);
zgnu(id3)=interp1(Tnodd,zgnu3,Tnod);

zgnl1=zgnl(id1);
zgnl2=zgnl(id2);
zgnl3=zgnl(id3);

zgnl(id1)=interp1(Tnodd,zgnl1,Tnod);
zgnl(id2)=interp1(Tnodd,zgnl2,Tnod);
zgnl(id3)=interp1(Tnodd,zgnl3,Tnod);

hmin=zeros(1,Nr+1);
yr=zeros(1,Nr+1);
fl=zeros(1,Nr+1);
zr=zeros(1,Nr+1);
alr=zeros(1,Nr+1);
z0=zeros(1,Nr+1);
for j=1:Nbe
    Fexte=zeros(16,1);
    ind=(Npe*(j-1)+1):(Npe*j+1);
    inde=(8*(j-1)+1):(8*(j+1));
    indg=(Nge*(j-1)+1):(Nge*j)+1;
    %======= calculation of groove and liner clearances =======%
    %nodal y and z displacements for element j
    Uey=Ug_nd([8*(j-1)+1:8*(j-1)+3,8*j+1:8*j+3]);
    Uez=Ug_nd([8*(j-1)+4:8*(j-1)+6,8*j+4:8*j+6]);
    %nodal twist for element j
    Uet=Ug_nd([8*j-1:8*j,8*j+7:8*j+8]);
    
    %ring radial displacement for element j
    yre=(Uey')*N*uref;
    %ring axial displacement for element j
    zre=(Uez')*N*href;
    zrg=(Uez')*N_g*href;
    zrgm1=zrg;
    %ring twist for element j
    alre=(Uet')*Nt*href/Rr;
    alrg=(Uet')*Ntg*href/Rr;
    alrgm1=alrg;
    %groove tilt for element j
    beta=beta_p*cos((Tnod(j)+Telt*x));
    betag=beta_p*cos((Tnod(j)+Telt*xg));
    betagm1=betag;
    
    %thermal tilting for element j
    beta_thu=tilt_thu*ones(size(x)); %upper flank
    beta_thug=tilt_thu*ones(size(xg)); %upper flank
    
    beta_thl=tilt_thl*ones(size(x)); %lower flank
    beta_thlg=tilt_thl*ones(size(xg)); %lower flank  
    %groove axial displacement due to tilt for element j
    dz_g=(Rr*cos((Tnod(j)+Telt*x))-off)*beta_p;
    dz_gg=(Rr*cos((Tnod(j)+Telt*xg))-off)*beta_p;
    dz_ggm1=dz_gg;
    
    %bore distortion
    id=3*(j-1)+1:3*(j+1);
    yb=(ybn(id))'*N*uref;
    
    %groove distortion
    zgu=(zgnu(id))'*N*href; %upper flank
    zgl=(zgnl(id))'*N*href; %lower flank
    zgug=(zgnu(id))'*N_g*href; %upper flank
    zglg=(zgnl(id))'*N_g*href; %lower flank
    
    templ=temp_liner(ind);
    %OCR clearance
    hrle=hrl(ind);
    Kpe=Kp(ind);
    Kfe=Kf(ind);
    
    hote=hot(indg);
    hobe=hob(indg);
    mu_oil=mu_oilv(indg);
    Pue=Pu(indg);
    Pde=Pd(indg);
    Pie=Pi(indg);
    isfuele=isfuel(ind);
    
    [fcgu,mcgu]=ring_groove_contact_vg(gcl,dz_g,zgu,zre,alre,beta,(beta_thu+thgt-thrt),y1u,y2u,sigmagt,z,Pk_rg,omega,1);
    [fcgl,mcgl]=ring_groove_contact_vg(gcl,dz_g,zgl,zre,alre,beta,(beta_thl-thgb+thrb),y1l,y2l,sigmagb,z,Pk_rg,omega,-1);
    [fcl,mcl,~]=ring_liner_contact_vg(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yre,alre,PR,sigmap,z,Pk_rl,omega);
    [hmine,~,fhl,ffl,mhl,z0l]=ring_liner_hydro_vg(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yre,alre,hrle,KOCR,Kpe,Kfe,Ph_OCR,ap,F0,b,c,d,e,rw,muUot,templ,isfuele,Vp,inputs.oil,sigmap,oilthreshold,viscosityfactor);
    [frgas,fzgas,mgas]=ring_gas_exg_vg(Pie,Pue,Pde,arm,aui,lfu,ali,lfl,hui,hli,huo,hlo,rbn,rb1,rb2,a11,a12,a21,a22,alrg);  
    [fgu,mgu]=ring_tgroove_og_vg(Pie,Pue,gcl,dz_gg,dz_ggm1,zgug,zrg,zrgm1,alrg,alrgm1,betag,betagm1,(beta_thug+thgt-thrt),y1u,y2u,hote,mu_oil,dt);
    [fgl,mgl]=ring_bgroove_og_vg(Pie,Pde,gcl,dz_gg,dz_ggm1,zglg,zrg,zrgm1,alrg,alrgm1,betag,betagm1,(beta_thlg-thgb+thrb),y1l,y2l,hobe,mu_oil,dt);   
    
    directionf=-sign(Vp);
    ffc=(-fcl)*fc*directionf;
    
    %z direction
    fi=Fi/2/pi/Rr;
    Ffi=sum(Le^4/Er/Iyy/href*N_int*fi*dx,2);
    
    fcz=ffc+fcgu+fcgl+ffl;
    fgz=fgu+fgl+fzgas;
    Fcz=Le^4/Er/Iyy/href*N_int*((fcz(1:end-1)+fcz(2:end))/2)'*dx;
    Fgz=Le^4/Er/Iyy/href*Ng_int*((fgz(1:end-1)+fgz(2:end))/2)'*dxg;
    
    %y direction
    fcy=-fcgu*thrt+fcgl*thrb+fcl+fhl;
    fgy=-fgu*thrt+fgl*thrb+frgas;
    
    Fcy=Le^4/Er/Izz/uref*N_int*((fcy(1:end-1)+fcy(2:end))/2)'*dx;  
    Fgy=Le^4/Er/Izz/uref*Ng_int*((fgy(1:end-1)+fgy(2:end))/2)'*dxg;
    
    %alpha direction
    mca=mcgu+fcgu*hui*thrt+mcgl+fcgl*hli*thrb+mcl+mhl+ffc*arm;
    mga=mgu+fgu*hui*thrt+mgl+fgl*hli*thrb+mgas;
    
    Mca=Rr*Le^2/Gr/Jt/href*Nt_int*((mca(1:end-1)+mca(2:end))/2)'*dx;
    Mga=Rr*Le^2/Gr/Jt/href*Ntg_int*((mga(1:end-1)+mga(2:end))/2)'*dxg;
    
    %coeff=alT*hr*ar^2*Rr^2*Telt^3/12/pi/Izz/uref;
    %Mgt=coeff*dTemp*dN_int*ones(Npe,1)*dx;
    %if j<floor((Nbe+1)/2)
    %   Mgt(1)=-Mgt(1);
    %   Mgt(2)=-Mgt(2);
    %   Mgt(3)=-Mgt(3);
    %end
    %if j+1<floor((Nbe+1)/2)
    %   Mgt(4)=-Mgt(4);
    %   Mgt(5)=-Mgt(5);
    %   Mgt(6)=-Mgt(6);
    %end
    Fexte([1:3,9:11])=Fcy+Fgy;
    %Fexte([1:3,9:11])=Fcy+Fgy+Mgt;
    Fexte([4:6,12:14])=Fcz+Fgz+Ffi;
    Fexte([7,8,15,16])=Mca+Mga;
    
    Fext_nd(inde)=Fext_nd(inde)+Fexte;
    
    hmin(ind)=hmine;
    z0(ind)=z0l;
    fl(ind)=fcl+fhl;
    zr(ind)=zre;
    yr(ind)=yre;
    alr(ind)=alre;
end
    MT_nd=Rr^3*Telt^3/Le/Er/Izz/uref*MT;
    %Fext_nd(2)=Fext_nd(2)-MT_nd;
    %Fext_nd(end-6)=Fext_nd(end-6)+MT_nd;
end