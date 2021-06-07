function [Fext_nd]=Ext_Forces_static_o(Ug_nd,inputs)
%Calculation of nodal resultants of external forces acting on the ring (static
%analysis)

%-------------- Retrieving parameters from input structure ---------------%
Db=inputs.bore.Db;
%Meshing
Telt=inputs.FEM.Telt; %element size (rad)
Tnod=inputs.FEM.Tnod; %node angular position (rad)
Nbe=inputs.FEM.Nbe; %number of elements
Nbnod=inputs.FEM.Nbnod; %number of nodes
Npe=inputs.FEM.Npe; %number of contact points within an element
uref=inputs.FEM.uref;
href=inputs.FEM.href;

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
thrt=inputs.ring.thrt;
thrb=inputs.ring.thrb;
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
sigmagu=inputs.piston.sigmag1t;
sigmagl=inputs.piston.sigmag1b;
Drg=inputs.piston.Drg;
hgi=inputs.piston.hgi;
gr_exp=inputs.piston.gr_exp;
thgt=inputs.piston.thgt;
thgb=inputs.piston.thgb;
hr=hui+thrt*aui+hli+thrb*ali;     %at centroid
hg=hgi+(Db/2-Drg/2-gr_exp-arm)*(thgt+thgb);   %at centroid

%liner
sigmap=inputs.liner.sigmap;
PR=inputs.liner.PR;

%contact
z=inputs.contact.z;
omega=inputs.contact.omega;
Pk_rg=inputs.contact.Pk_rg;
Pk_rl=inputs.contact.Pk_rl;

%secondary motion
beta_p=0; %piston tilt - current step
off=0; %pin offset
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

%Shape functions for twist displacements
Nt1=1-3*x.^2+2*x.^3;
Nt2=x-2*x.^2+x.^3;
Nt3=3*x.^2-2*x.^3;
Nt4=-x.^2+x.^3;

Nt=[Nt1;Nt2;Nt3;Nt4];
Nt_int=(Nt(:,1:end-1)+Nt(:,2:end))/2;

%Initializing vector of external forces
Fext_nd=zeros(8*Nbnod,1);

%Bore distortion at node loactions
[ybn]=bore_dist_static_o(inputs);

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
[zgnu,zgnl]=groove_dist_static_o(inputs);

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

for j=1:Nbe
    Fexte=zeros(16,1);
    indg=(8*(j-1)+1):(8*(j+1));
    %======= calculation of groove and liner clearances =======%
    %nodal y and z displacements for element j
    Uey=Ug_nd([8*(j-1)+1:8*(j-1)+3,8*j+1:8*j+3]);
    Uez=Ug_nd([8*(j-1)+4:8*(j-1)+6,8*j+4:8*j+6]);
    %nodal twist for element j
    Uet=Ug_nd([8*j-1:8*j,8*j+7:8*j+8]);
    
    %ring radial displacement for element j
    yr=(Uey')*N*uref;
    %ring axial displacement for element j
    zr=(Uez')*N*href;
    %ring twist for element j
    alr=(Uet')*Nt*href/Rr;
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
    
    [fcgu,mcgu]=ring_groove_contact_o(gcl,dz_g,zgu,zr,alr,beta,(beta_thu+thgt-thrt),y1u,y2u,sigmagu,z,Pk_rg,omega,1);
    [fcgl,mcgl]=ring_groove_contact_o(gcl,dz_g,zgl,zr,alr,beta,(beta_thl-thgb+thrb),y1l,y2l,sigmagl,z,Pk_rg,omega,-1);
    
    [fcl,mcl,~]=ring_liner_contact_o(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yr,alr,PR,sigmap,z,Pk_rl,omega);
    
    %z direction
    Fcgu=Le^4/Er/Iyy/href*N_int*((fcgu(1:end-1)+fcgu(2:end))/2)'*dx;
    Fcgl=Le^4/Er/Iyy/href*N_int*((fcgl(1:end-1)+fcgl(2:end))/2)'*dx;
   
    %y direction
    Fcl=Le^4/Er/Izz/uref*N_int*((fcl(1:end-1)+fcl(2:end))/2)'*dx;  %in no dimension
    Fcgur=Le^4/Er/Izz/uref*N_int*(-(fcgu(1:end-1)+fcgu(2:end))*thrt/2)'*dx;
    Fcglr=Le^4/Er/Izz/uref*N_int*((fcgl(1:end-1)+fcgl(2:end))*thrb/2)'*dx;
    
    %alpha direction
    Mcgu=Rr*Le^2/Gr/Jt/href*Nt_int*((mcgu(1:end-1)+mcgu(2:end))/2+...
        (fcgu(1:end-1)+fcgu(2:end))/2*hui*thrt)'*dx;
    Mcgl=Rr*Le^2/Gr/Jt/href*Nt_int*((mcgl(1:end-1)+mcgl(2:end))/2+...
        (fcgl(1:end-1)+fcgl(2:end))/2*hli*thrb)'*dx;
    Mcl=Rr*Le^2/Gr/Jt/href*Nt_int*((mcl(1:end-1)+mcl(2:end))/2)'*dx;
    
    Fexte([1:3,9:11])=Fcl+Fcgur+Fcglr;
    Fexte([4:6,12:14])=Fcgu+Fcgl;
    Fexte([7,8,15,16])=Mcgu+Mcgl+Mcl;
    
    Fext_nd(indg)=Fext_nd(indg)+Fexte;
end
    Fext_nd(2)=Fext_nd(2);
    Fext_nd(end-6)=Fext_nd(end-6);
end