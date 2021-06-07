function [Fext_nd,yr,zr,alr,hmin,fl,z0]=Ext_Forces_static_n4g(Ug_nd,inputs)
%Calculation of nodal resultants of external forces acting on the ring (static
%analysis)

%-------------- Retrieving parameters from input structure ---------------%
Db=inputs.bore.Db;

%engine
Vp=inputs.engine.Vp;

%oil data
% zk=inputs.oil.zk;
% temp1=inputs.oil.temp1;
% temp2=inputs.oil.temp2;
% rho_oil=inputs.oil.rho_oil;

%Meshing
Telt=inputs.FEM.Telt; %element size (rad)
Tnod=inputs.FEM.Tnod; %node angular position (rad)
Nbe=inputs.FEM.Nbe; %number of elements
Nbnod=inputs.FEM.Nbnod; %number of nodes
Nr=inputs.FEM.Nr;
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
arm=inputs.ring.arm;
ali=inputs.ring.ali;
alo=inputs.ring.alo;
aui=inputs.ring.aui;
ar=arm+max(ali+aui);
rb1=inputs.ring.rb1;
rb2=inputs.ring.rb2;
rbn=inputs.ring.rbn;
a10=inputs.ring.a10;
a11=inputs.ring.a11;
a12=inputs.ring.a12;
a20=inputs.ring.a20;
a21=inputs.ring.a21;
a22=inputs.ring.a22;
hui=inputs.ring.hui;
hli=inputs.ring.hli;
hlo=inputs.ring.hlo;
thrt=inputs.ring.thrt;
thrb=inputs.ring.thrb;
alp=inputs.ring.alp;
isi=inputs.ring.isi;
Tg=inputs.ring.Tg;

%groove
gcl=inputs.groove.gcl;
y1u=inputs.groove.y1u;
y2u=inputs.groove.y2u;
y1l=inputs.groove.y1l;
y2l=inputs.groove.y2l;
tilt_thu=inputs.groove.tilt_thu; %groove thermal tilting - upper flank
tilt_thl=inputs.groove.tilt_thl; %groove thermal tilting - lower flank

%piston
sigmagt=inputs.piston.sigmag1t;
sigmagb=inputs.piston.sigmag1b;
%Drg=inputs.piston.Drg;
hgi=inputs.piston.hgi;
%gr_exp=inputs.piston.gr_exp;
thgt=inputs.piston.thgt;
thgb=inputs.piston.thgb;
hr=hui+thrt*aui+hli+thrb*ali;     %at centroid
%hg=hgi+(Db/2-Drg/2-gr_exp-arm)*(thgt+thgb);   %at centroid

%liner
sigmap=inputs.liner.sigmap;
PR=inputs.liner.PR;

%contact
z=inputs.contact.z;
omega=inputs.contact.omega;
Pk_rg=inputs.contact.Pk_rg;
Pk_rl=inputs.contact.Pk_rl;
fc=inputs.contact.fc_dry;

%secondary motion
beta_p=inputs.pistonSM.beta_p; %piston tilt - current step
off=inputs.pistonSM.off; %pin offset
%-------------------------------------------------------------------------%

%discretization inside the element
dx=1/Npe;
x=0:dx:1;

%Shape functions for z displacements
N1=1-10*x.^3+15*x.^4-6*x.^5;
N2=x-6*x.^3+8*x.^4-3*x.^5;
N3=(x.^2)/2-3*(x.^3)/2+3*(x.^4)/2-(x.^5)/2;
N4=10*x.^3-15*x.^4+6*x.^5;
N5=-4*x.^3+7*x.^4-3*x.^5;
N6=(x.^3)/2-x.^4+x.^5/2;
N=[N1;N2;N3;N4;N5;N6];
N_int=(N(:,1:end-1)+N(:,2:end))/2;

%Initializing vector of external forces
Fext_nd=zeros(12*Nbnod,1);

%Bore distortion at node loactions
[ybn]=bore_dist_static_n4g(inputs);

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
[zgnu,zgnl]=groove_dist_static_n4g(inputs);

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
    Fexte=zeros(24,1);
    ind=(Npe*(j-1)+1):(Npe*j+1);
    inde=(12*(j-1)+1):(12*(j+1));
    %indg=(Nge*(j-1)+1):(Nge*j)+1;
    %======= calculation of groove and liner clearances =======%
    %nodal y and z displacements for element j
    Uey=Ug_nd([12*(j-1)+1:12*(j-1)+3,12*j+1:12*j+3]);
    Uez=Ug_nd([12*(j-1)+4:12*(j-1)+6,12*j+4:12*j+6]);
    %nodal twist for element j
    Uet=Ug_nd([12*(j-1)+7:12*(j-1)+9,12*j+7:12*j+9]);
    
    %ring radial displacement for element j
    yre=(Uey')*N*uref;
    %ring axial displacement for element j
    zre=(Uez')*N*href;
    
    %ring twist for element j
    alre=(Uet')*N*href/Rr;
    
    %groove tilt for element j
    beta=beta_p*cos((Tnod(j)+Telt*x));
   
    %thermal tilting for element j
    beta_thu=tilt_thu*ones(size(x)); %upper flank
   
    beta_thl=tilt_thl*ones(size(x)); %lower flank
    
    %groove axial displacement due to tilt for element j
    dz_g=(Rr*cos((Tnod(j)+Telt*x))-off)*beta_p;
    
    %bore distortion
    id=3*(j-1)+1:3*(j+1);
    yb=(ybn(id))'*N*uref;
    
    %groove distortion
    zgu=(zgnu(id))'*N*href; %upper flank
    zgl=(zgnl(id))'*N*href; %lower flank
    
    [fcgu,mcgu]=ring_groove_contact_n4g(gcl,dz_g,zgu,zre,alre,alp,beta,(beta_thu+thgt-thrt),y1u,y2u,sigmagt,z,Pk_rg,omega,1,thrt,thrb);
    [fcgl,mcgl]=ring_groove_contact_n4g(gcl,dz_g,zgl,zre,alre,alp,beta,(beta_thl-thgb+thrb),y1l,y2l,sigmagb,z,Pk_rg,omega,-1,thrt,thrb);
    [fcl,mcl,hmine,z01]=ring_liner_contact_n4g(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yre,alre,alp,PR,sigmap,z,Pk_rl,omega);
  
    fcz=fcgu+fcgl;
    Fcz=Le^4/Er/Iyy/href*N_int*((fcz(1:end-1)+fcz(2:end))/2)'*dx;
    
    %y direction
    fcy=-fcgu*thrt+fcgl*thrb+fcl;

    Fcy=Le^4/Er/Izz/uref*N_int*((fcy(1:end-1)+fcy(2:end))/2)'*dx;  
    
    %alpha direction
    mca=mcgu+fcgu*hui*thrt+mcgl+fcgl*hli*thrb+mcl;
    
    Mca=Rr*Le^2/Gr/Jt/href*N_int*((mca(1:end-1)+mca(2:end))/2)'*dx;
    
    Fexte([1:3,13:15])=Fcy;
    Fexte([4:6,16:18])=Fcz;
    Fexte([7:9,19:21])=Mca;
    Fexte([10,22])=(-gcl/2+omega*sigmagb)/href*ones(2,1);
    
    Fext_nd(inde)=Fext_nd(inde)+Fexte;
    
    hmin(ind)=hmine;
    z0(ind)=z01;
    fl(ind)=fcl;
    zr(ind)=zre;
    yr(ind)=yre;
    alr(ind)=alre;
end
    %MT_nd=Rr^3*Telt^3/Le/Er/Izz/uref*MT;
    %Fext_nd(2)=Fext_nd(2)-MT_nd;
    %Fext_nd(end-6)=Fext_nd(end-6)+MT_nd;
end