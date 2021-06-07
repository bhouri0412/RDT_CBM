function [Ja_stac]=Jacobian_Estac_vg(Ug_nd,inputs)
%Calculates the Jacobian of the error function of the ring dynamics

%-------------- Retrieving parameters from input structure ---------------%
Db=inputs.bore.Db;
temp=inputs.bore.templ;
dt=1;

%engine
Vp=inputs.engine.Vp;
Pu=inputs.engine.Pu;
Pd=inputs.engine.Pd;
Pi=inputs.engine.Pi;
%Meshing
Telt=inputs.FEM.Telt; %element size (rad)
Tnod=inputs.FEM.Tnod; %node angular position (rad)
Nbe=inputs.FEM.Nbe; %number of elements
Nbnod=inputs.FEM.Nbnod; %number of nodes
Npe=inputs.FEM.Npe; %number of contact points within an element
Nr=inputs.FEM.Nr;
Ng=inputs.FEM.Ng;
uref=inputs.FEM.uref;
href=inputs.FEM.href;
Kg_nd=inputs.FEM.Kg_nd;
Nge=inputs.FEM.Nge;

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
%dTemp=inputs.ring.dTemp;
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
lfu=inputs.groove.lfu;
lfl=inputs.groove.lfl;
tilt_thu=inputs.groove.tilt_thu; %groove thermal tilting - upper flank
tilt_thl=inputs.groove.tilt_thl; %groove thermal tilting - lower flank

mu_oilv=inputs.oil.mu_oil;

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
%hot=4e-6;
%hob=4e-6;

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
%secondary motion
beta_p=inputs.pistonSM.beta_p; %piston tilt - current step
off=inputs.pistonSM.off; %pin offset
%-------------------------------------------------------------------------%
shiftstep=round(Tg/2/pi*Nr);
temp_liner=temp;

%shift to ring reference
isfuel=isfuela;
hrls=hrls.*(hrls>=(minoil))+minoil*(hrls<(minoil));
hrl=2*hrls;
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
Ja_F=zeros(8*Nbnod,8*Nbnod);

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

for j=1:Nbe
    JaFe=zeros(16,16);
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
    yr=(Uey')*N*uref;
    %ring axial displacement for element j
    zr=(Uez')*N*href;
    zrg=(Uez')*N_g*href;
    zrgm1=zrg;
    %ring twist for element j
    alr=(Uet')*Nt*href/Rr;
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
    
    [dfcldy,dfcldalr,dmcldy,dmcldalr]=diff_ring_liner_contact_vg(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yr,alr,PR,sigmap,z,Pk_rl,omega);
    [dfhldy,dffldy,dmhldy,dfhldalr,dffldalr,dmhldalr]=diff_ring_liner_hydro_vg(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yr,alr,hrle,KOCR,Kpe,Kfe,Ph_OCR,ap,F0,b,c,d,e,rw,muUot,templ,isfuele,Vp,inputs.oil,sigmap,oilthreshold,viscosityfactor);
    [dfcgudz,dfcgudalr,dmcgudz,dmcgudalr]=diff_ring_groove_contact_vg(gcl,dz_g,zgu,zr,alr,beta,(beta_thu+thgt-thrt),y1u,y2u,sigmagu,z,Pk_rg,omega,1);
    [dfcgldz,dfcgldalr,dmcgldz,dmcgldalr]=diff_ring_groove_contact_vg(gcl,dz_g,zgl,zr,alr,beta,(beta_thl-thgb+thrb),y1l,y2l,sigmagl,z,Pk_rg,omega,-1);
    [dfgudz,dfgudalr,dmgudz,dmgudalr,~,~,~,~]=diff_ring_tgroove_og_vg(Pie,Pue,gcl,dz_gg,dz_ggm1,zgug,zrg,zrgm1,alrg,alrgm1,betag,betagm1,(beta_thug+thgt-thrt),y1u,y2u,hote,mu_oil,dt);
    [dfgldz,dfgldalr,dmgldz,dmgldalr,~,~,~,~]=diff_ring_bgroove_og_vg(Pie,Pde,gcl,dz_gg,dz_ggm1,zglg,zrg,zrgm1,alrg,alrgm1,betag,betagm1,(beta_thlg-thgb+thrb),y1l,y2l,hobe,mu_oil,dt);
    [dfrgasdalr,dmgasdalr]=diff_ring_gas_exg_vg(Pue,Pde,rbn,rb1,rb2,a11,a12,a21,a22,alrg); 
    
    %r direction
    dfcydy=dfcldy+dfhldy;
    dfcydz=-dfcgudz*thrt+dfcgldz*thrb;
    dfgydz=-dfgudz*thrt+dfgldz*thrb;
    dfcydalr=dfcldalr-dfcgudalr*thrt+dfcgldalr*thrb+dfhldalr;
    dfgydalr=-dfgudalr*thrt+dfgldalr*thrb+dfrgasdalr;
    
    dFcydy=Le^4/Er/Izz*(repmat((dfcydy(1:end-1)+dfcydy(2:end))/2,6,1).*N_int*dx)*N_int';
    
    dFcydz=Le^4*href/uref/Er/Izz*(repmat(-(dfcydz(1:end-1)+dfcydz(2:end))/2,6,1).*N_int*dx)*N_int';
    dFgydz=Le^4*href/uref/Er/Izz*(repmat((dfgydz(1:end-1)+dfgydz(2:end))/2,6,1).*Ng_int*dxg)*Ng_int';
    
    dFcydalr=Le^4*href/Er/Izz/uref/Rr*(repmat((dfcydalr(1:end-1)+dfcydalr(2:end))/2,6,1).*N_int*dx)*Nt_int';
    dFgydalr=Le^4*href/Er/Izz/uref/Rr*(repmat(-(dfgydalr(1:end-1)+dfgydalr(2:end))/2,6,1).*Ng_int*dxg)*Ntg_int';
    
    %z direction     
    directionf=-sign(Vp);
    dffcdy=-fc*directionf*dfcldy;
    dffcdalr=-fc*directionf*dfcldalr;
    
    dfczdy=dffcdy+dffldy;
    
    dfczdz=dfcgudz+dfcgldz;
    dfgzdz=dfgudz+dfgldz;
    
    dfczdalr=dffcdalr+dffldalr+dfcgudalr+dfcgldalr;
    dfgzdalr=dfgudalr+dfgldalr;
    
    dFffdy=Le^4*uref/Er/Iyy/href*(repmat((dfczdy(1:end-1)+dfczdy(2:end))/2,6,1).*N_int*dx)*N_int';
    
    dFczdz=Le^4/Er/Iyy*(repmat((dfczdz(1:end-1)+dfczdz(2:end))/2,6,1).*N_int*dx)*N_int';
    dFgzdz=Le^4/Er/Iyy*(repmat((dfgzdz(1:end-1)+dfgzdz(2:end))/2,6,1).*Ng_int*dxg)*Ng_int';
    
    dFczdalr=Le^4/Er/Iyy/Rr*(repmat((dfczdalr(1:end-1)+dfczdalr(2:end))/2,6,1).*N_int*dx)*Nt_int';
    dFgzdalr=Le^4/Er/Iyy/Rr*(repmat((dfgzdalr(1:end-1)+dfgzdalr(2:end))/2,6,1).*Ng_int*dxg)*Ntg_int';
    
    %alpha direction
    dmcdy=dmcldy+dmhldy+dffcdy*arm;
    dmcdz=dmcgudz+dfcgudz*hui*thrt+dmcgldz+dfcgldz*hli*thrb;
    dmgdz=dmgudz+dfgudz*hui*thrt+dmgldz+dfgldz*hli*thrb;
    dmcdalr=dmcldalr+dmhldalr+dffcdalr*arm+dmcgudalr+dfcgudalr*hui*thrt+dmcgldalr+dfcgldalr*hli*thrb;
    dmgdalr=dmgudalr+dfgudalr*hui*thrt+dmgldalr+dfgldalr*hli*thrb+dmgasdalr;
    
    dMcdy=Rr*Le^2*uref/Gr/Jt/href*(repmat((dmcdy(1:end-1)+dmcdy(2:end))/2,4,1).*Nt_int*dx)*N_int';
    
    dMcdz=Rr*Le^2/Gr/Jt*(repmat(((dmcdz(1:end-1)+dmcdz(2:end))/2),4,1).*Nt_int*dx)*N_int';
    dMgdz=Rr*Le^2/Gr/Jt*(repmat(((dmgdz(1:end-1)+dmgdz(2:end))/2),4,1).*Ntg_int*dxg)*Ng_int';
    
    dMcdalr=Le^2/Gr/Jt*(repmat(((dmcdalr(1:end-1)+dmcdalr(2:end))/2),4,1).*Nt_int*dx)*Nt_int'; 
    dMgdalr=Le^2/Gr/Jt*(repmat(((dmgdalr(1:end-1)+dmgdalr(2:end))/2),4,1).*Ntg_int*dxg)*Ntg_int';   
    %dMgldalr=Le^2/Gr/Jt*(repmat(((dmcgldalr(1:end-1)+dmcgldalr(2:end))/2+(dfcgldalr(1:end-1)+dfcgldalr(2:end))*hli*thrb/2),4,1).*Nt_int*dx)*Nt_int'; 
    
    JaFe([1:3,9:11],[1:3,9:11])=dFcydy;
    JaFe([1:3,9:11],[4:6,12:14])=dFcydz+dFgydz;
    JaFe([1:3,9:11],[7:8,15:16])=dFcydalr+dFgydalr;
    
    JaFe([4:6,12:14],[1:3,9:11])=dFffdy;
    JaFe([4:6,12:14],[4:6,12:14])=dFgzdz+dFczdz;
    JaFe([4:6,12:14],[7:8,15:16])=dFgzdalr+dFczdalr;
    
    JaFe([7:8,15:16],[1:3,9:11])=dMcdy;
    JaFe([7:8,15:16],[4:6,12:14])=dMcdz+dMgdz;
    JaFe([7:8,15:16],[7:8,15:16])=dMcdalr+dMgdalr;
    
    Ja_F(inde,inde)=Ja_F(inde,inde)+JaFe;  
    
end

Ja_stac=Kg_nd-Ja_F+1e-12;

end