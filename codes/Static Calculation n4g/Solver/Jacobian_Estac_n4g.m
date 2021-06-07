function [Ja_stac]=Jacobian_Estac_n4g(Ug_nd,inputs)
%Calculates the Jacobian of the error function of the ring dynamics

%-------------- Retrieving parameters from input structure ---------------%
Db=inputs.bore.Db;

%engine
Vp=inputs.engine.Vp;

%Meshing
Telt=inputs.FEM.Telt; %element size (rad)
Tnod=inputs.FEM.Tnod; %node angular position (rad)
Nbe=inputs.FEM.Nbe; %number of elements
Nbnod=inputs.FEM.Nbnod; %number of nodes
Npe=inputs.FEM.Npe; %number of contact points within an element
uref=inputs.FEM.uref;
href=inputs.FEM.href;
Kg_nd=inputs.FEM.Kg_nd;

%ring
Er=inputs.ring.Er;
%alT=inputs.ring.alT;
Gr=inputs.ring.Gr;
Izz=inputs.ring.Izz;
Iyy=inputs.ring.Iyy;
Jt=inputs.ring.Jt;
Rr=inputs.ring.Rr; %ring neutral fiber radius
Le=inputs.ring.Le;
%arm=inputs.ring.arm;
% ali=inputs.ring.ali;
% aui=inputs.ring.aui;
%ar=arm+max(ali+aui);
% huo=inputs.ring.huo;
% hlo=inputs.ring.hlo;
arm=inputs.ring.arm; %ring cross section centroid coordinate
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
alp=inputs.ring.alp;
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
% Drg=inputs.piston.Drg;
% hgi=inputs.piston.hgi;
% gr_exp=inputs.piston.gr_exp;
thgt=inputs.piston.thgt;
thgb=inputs.piston.thgb;
% hr=hui+thrt*aui+hli+thrb*ali;     %at centroid
% hg=hgi+(Db/2-Drg/2-gr_exp-arm)*(thgt+thgb);   %at centroid

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
Ja_F=zeros(12*Nbnod,12*Nbnod);

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

for j=1:Nbe
    JaFe=zeros(24,24);
    ind=(Npe*(j-1)+1):(Npe*j+1);
    inde=(12*(j-1)+1):(12*(j+1));
    %======= calculation of groove and liner clearances =======%
    %nodal y and z displacements for element j
    Uey=Ug_nd([12*(j-1)+1:12*(j-1)+3,12*j+1:12*j+3]);
    Uez=Ug_nd([12*(j-1)+4:12*(j-1)+6,12*j+4:12*j+6]);
    %nodal twist for element j
    Uet=Ug_nd([12*(j-1)+7:12*(j-1)+9,12*j+7:12*j+9]);
    
    %ring radial displacement for element j
    yr=(Uey')*N*uref;
    %ring axial displacement for element j
    zr=(Uez')*N*href;
    
    %ring twist for element j
    alr=(Uet')*N*href/Rr;
    
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
    
    [dfcldy,dfcldalr,dmcldy,dmcldalr]=diff_ring_liner_contact_n4g(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yr,alr,alp,PR,sigmap,z,Pk_rl,omega);
    [dfcgudz,dfcgudalr,dmcgudz,dmcgudalr]=diff_ring_groove_contact_n4g(gcl,dz_g,zgu,zr,alr,alp,beta,(beta_thu+thgt-thrt),y1u,y2u,sigmagu,z,Pk_rg,omega,1,thrt,thrb);
    [dfcgldz,dfcgldalr,dmcgldz,dmcgldalr]=diff_ring_groove_contact_n4g(gcl,dz_g,zgl,zr,alr,alp,beta,(beta_thl-thgb+thrb),y1l,y2l,sigmagl,z,Pk_rg,omega,-1,thrt,thrb);
    
    %r direction
    dfcydy=dfcldy;
    dfcydz=-dfcgudz*thrt+dfcgldz*thrb;
    dfcydalr=dfcldalr-dfcgudalr*thrt+dfcgldalr*thrb;
    
    dFcydy=Le^4/Er/Izz*(repmat((dfcydy(1:end-1)+dfcydy(2:end))/2,6,1).*N_int*dx)*N_int';
    
    dFcydz=Le^4*href/uref/Er/Izz*(repmat(-(dfcydz(1:end-1)+dfcydz(2:end))/2,6,1).*N_int*dx)*N_int';
    
    dFcydalr=Le^4*href/Er/Izz/uref/Rr*(repmat((dfcydalr(1:end-1)+dfcydalr(2:end))/2,6,1).*N_int*dx)*N_int';
   
    %z direction     
   
    dfczdz=dfcgudz+dfcgldz;
    
    dfczdalr=dfcgudalr+dfcgldalr;
        
    dFczdz=Le^4/Er/Iyy*(repmat((dfczdz(1:end-1)+dfczdz(2:end))/2,6,1).*N_int*dx)*N_int';
   
    dFczdalr=Le^4/Er/Iyy/Rr*(repmat((dfczdalr(1:end-1)+dfczdalr(2:end))/2,6,1).*N_int*dx)*N_int';
    
    %alpha direction
    dmcdy=dmcldy;
    dmcdz=dmcgudz+dfcgudz*hui*thrt+dmcgldz+dfcgldz*hli*thrb;
    dmcdalr=dmcldalr+dmcgudalr+dfcgudalr*hui*thrt+dmcgldalr+dfcgldalr*hli*thrb;
   
    dMcdy=Rr*Le^2*uref/Gr/Jt/href*(repmat((dmcdy(1:end-1)+dmcdy(2:end))/2,6,1).*N_int*dx)*N_int';
    
    dMcdz=Rr*Le^2/Gr/Jt*(repmat(((dmcdz(1:end-1)+dmcdz(2:end))/2),6,1).*N_int*dx)*N_int';

    dMcdalr=Le^2/Gr/Jt*(repmat(((dmcdalr(1:end-1)+dmcdalr(2:end))/2),6,1).*N_int*dx)*N_int'; 
    
    JaFe([1:3,13:15],[1:3,13:15])=dFcydy;
    JaFe([1:3,13:15],[4:6,16:18])=dFcydz;
    JaFe([1:3,13:15],[7:9,19:21])=dFcydalr;
    
    JaFe([4:6,16:18],[4:6,16:18])=dFczdz;
    JaFe([4:6,16:18],[7:9,19:21])=dFczdalr;
    
    JaFe([7:9,19:21],[1:3,13:15])=dMcdy;
    JaFe([7:9,19:21],[4:6,16:18])=dMcdz;
    JaFe([7:9,19:21],[7:9,19:21])=dMcdalr;
    
    Ja_F(inde,inde)=Ja_F(inde,inde)+JaFe;  
    
end

Ja_stac=Kg_nd-Ja_F+1e-12;

end