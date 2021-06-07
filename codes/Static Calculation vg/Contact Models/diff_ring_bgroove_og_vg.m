function [dfgldz,dfgldalr,dmgldz,dmgldalr,dfgldpd,dfgldpi,dmgldpd,dmgldpi]...
    =diff_ring_bgroove_og_vg(Pi,Pd,gcl,dz_g,dz_go,zg,zr,zro,alr,alro,beta,betao,beta_th,y1,y2,hoil,mu_oil,dt)

Npe=length(zr);
Nby=50;

dfgldz=zeros(1,Npe);
dfgldalr=zeros(1,Npe);
dmgldz=zeros(1,Npe);
dmgldalr=zeros(1,Npe);
dfgldpd=zeros(1,Npe);
dfgldpi=zeros(1,Npe);
dmgldpd=zeros(1,Npe);
dmgldpi=zeros(1,Npe);

h0=gcl/2-zg-dz_g+zr;
al=(alr-beta-beta_th);
% al=al.*(abs(al)>1e-12)+1e-12*sign(al).*(abs(al)<1e-12);
% al=al.*(al~=0)+1e-12*(al==0);

h0o=gcl/2-zg-dz_go+zro;
alo=(alro-betao-beta_th);
% alo=alo.*(abs(alo)>1e-12)+1e-12*sign(alo).*(abs(alo)<1e-12);
% alo=alo.*(alo~=0)+1e-12*(alo==0);

minh=(h0+y2*al).*(al<0) ...
    +(h0+y1*al).*(al>0);

ifgas=minh>hoil;
ngas=nnz(ifgas);
%gas region
if ngas
    ifgas1=ifgas&(abs(al)>1e-8);
    ifgas2=ifgas&(abs(al)<=1e-8);
    ngas1=nnz(ifgas1);
    ngas2=nnz(ifgas2);
    
    if ngas1
        h0g1=h0(ifgas1)-hoil(ifgas1);
        alg=al(ifgas1);
        hc1=h0g1+alg*y1;
        hc1=hc1.*(hc1>=1e-12)+1e-12*(hc1<1e-12);
        hc2=h0g1+alg*y2;
        hc2=hc2.*(hc2>=1e-12)+1e-12*(hc2<1e-12);
        Pig1=Pi(ifgas1);
        Pdg1=Pd(ifgas1);
        
        h12p=(hc1+hc2);
        h12d=(hc1-hc2);
        h12p2=h12p.^2;
        h12d2=h12d.^2;
        
        fcom=h12d./alg./h12p;
        mcom=2*(hc1.*hc2).^2.*log(hc1./hc2);
        mcomd=1./(alg.^2.*h12p.*h12d);
        dfcom=h12d2./h12p2;
        dmcom=(h12d.*(h0g1.*h12d2+3*hc1.*hc2.*h12p)-2*hc1.*hc2.*(hc1.^2+hc1.*hc2+hc2.^2).*log(hc1./hc2))...
            ./(h12d.*h12p2);
        
        dfgldz(ifgas1)=(Pig1-Pdg1).*dfcom./alg;
        dfgldalr(ifgas1)=(-Pig1+Pdg1).*h0g1.*dfcom./(alg.^2);
        
        dmgldz(ifgas1)=(-Pig1+Pdg1).*dmcom./(alg.^2);
        dmgldalr(ifgas1)=(Pig1-Pdg1).*h0g1.*dmcom./(alg.^3);
        
        dfgldpi(ifgas1)=-hc1.*fcom;
        dfgldpd(ifgas1)=-hc2.*fcom;
        dmgldpi(ifgas1)=(mcom+hc1.*h12d.*(2*h0g1.*h12d-hc1.*h12p)).*mcomd/2;
        dmgldpd(ifgas1)=(-mcom+hc2.*h12d.*(2*h0g1.*h12d+hc2.*h12p)).*mcomd/2;
    end
    if ngas2
        h0g2=h0(ifgas2)-hoil(ifgas2);
        Pig2=Pi(ifgas2);
        Pdg2=Pd(ifgas2);
        comn=(Pdg2-Pig2).*(y1-y2)^2./h0g2/4;
        %dfgldz(ifgas2)=0;
        dfgldalr(ifgas2)=comn;
        %dmgldz(ifgas2)=0;
        dmgldalr(ifgas2)=comn*(y1+y2);
        
        dfgldpi(ifgas2)=1/2*(y2-y1);
        dfgldpd(ifgas2)=1/2*(y2-y1);
        dmgldpi(ifgas2)=1/6*(y2-y1)*(2*y1+y2);
        dmgldpd(ifgas2)=1/6*(y2-y1)*(y1+2*y2);
    end
end

%oil region%
ifoil=minh<hoil;
noil=nnz(ifoil);
if noil
    h0oil=h0(ifoil);
    h0ooil=h0o(ifoil);
    
    aloil=al(ifoil);
    aloil=aloil.*(abs(aloil)>1e-8)+1e-8*sign(aloil).*(abs(aloil)<1e-8);
    aloil=aloil.*(aloil~=0)+1e-8*(aloil==0);
    
    alooil=alo(ifoil);
    alooil=alooil.*(abs(alooil)>1e-8)+1e-8*sign(alooil).*(abs(alooil)<1e-8);
    alooil=alooil.*(alooil~=0)+1e-8*(alooil==0);
    
    Pdoil=Pd(ifoil);
    Pioil=Pi(ifoil);
    mu_oilo=mu_oil(ifoil);
    
    yo=(hoil(ifoil)-h0oil)./aloil;
    yo=min(yo,y2-1e-6);
    yo=max(yo,y1+1e-6);
    yb=y2*(aloil<0)+yo.*(aloil>0);
    yb=((yb-y1)>1e-8).*yb+(y1+1e-8).*((yb-y1)<1e-8);    %yb cannot smaller than y1
    ya=y1*(aloil>0)+yo.*(aloil<0);    
    ya=((y2-ya)>1e-8).*ya+(y2-1e-8).*((y2-ya)<1e-8);  %ya cannot excess y2
    ya2=repmat(ya,Nby+1,1);
    
    dyoil=(yb-ya)/Nby;
    dyoil2=repmat(dyoil,Nby,1);
    yyoil=repmat(ya,Nby+1,1)+repmat([0:Nby]',1,noil).*repmat(dyoil,Nby+1,1);
    h0_2oil=repmat(h0oil,Nby+1,1);
    al_2oil=repmat(al(ifoil),Nby+1,1);
    hcl=h0_2oil+al_2oil.*yyoil;
    hcl=hcl.*(hcl>1e-7)+1e-7*(hcl<1e-7);
    %hclcheck=hcl>1e-8;
    hcl3=hcl.^(-3);
    hcl4m3=3*hcl.^(-4);
    yoad=yyoil-ya2;
    yoad2=yyoil.^2-ya2.^2;
    
    h0_2oilo=repmat(h0ooil,Nby+1,1);
    al_2oilo=repmat(alooil,Nby+1,1);
    hclo=h0_2oilo+al_2oilo.*yyoil;
    hclo=hclo.*(hclo>1e-7)+1e-7*(hclo<1e-7);
    
    dhdt=(hcl-hclo)/dt;
    dhdt=dhdt.*(dhdt<0);
    
    h2=h0oil+aloil*y2;
    h2=h2.*(h2>1e-7)+1e-7*(h2<1e-7);
    h1=h0oil+aloil*y1;
    h1=h1.*(h1>1e-7)+1e-7*(h1<1e-7);
    hb=hcl(end,:);
    ha=hcl(1,:);
    %hao=hclo(1,:);
    hb3=hb.^(-3);
    ha3=ha.^(-3);
    
    dybdz=-1./(aloil).*(aloil>0).*(h2>hoil(ifoil));
    dybda=-(hoil(ifoil)-h0oil)./(aloil.^2).*(aloil>0).*(h2>hoil(ifoil));
    dyadz=-1./(aloil).*(aloil<0).*(h1>hoil(ifoil));
    dyada=-(hoil(ifoil)-h0oil)./(aloil.^2).*(aloil<0).*(h1>hoil(ifoil));
    dyadz2=repmat(dyadz,Nby+1,1);
    dyada2=repmat(dyada,Nby+1,1);
    dhat2=repmat(dhdt(1,:),Nby+1,1);
    
    J=cumsum((dhdt(1:end-1,:)+dhdt(2:end,:))/2).*dyoil2;
    J=[zeros(1,noil);J];
    dt2=1/dt*ones(Nby+1,noil).*(dhdt<0);
    ddJdz=cumsum((dt2(1:end-1,:)+dt2(2:end,:))/2).*dyoil2;
    ddJdz=[zeros(1,nnz(ifoil));ddJdz];
    dJdz=-dyadz2.*dhat2+ddJdz;
    dty2=yyoil/dt.*(dhdt<0);
    ddJda=cumsum((dty2(1:end-1,:)+dty2(2:end,:))/2).*dyoil2;
    ddJda=[zeros(1,nnz(ifoil));ddJda];
    dJda=-dyada2.*dhat2+ddJda;
    
    Jb=J(end,:);
    
    dJ0=hcl3;
    ddJ0dz=-hcl4m3;
    ddJ0da=-yyoil.*hcl4m3;
    
    ddJJ1dz=(J.*(-dyadz2)+yoad.*dJdz).*hcl3...
        -yoad.*J.*hcl4m3;
    ddJJ1da=(-dyada2.*J+yoad.*dJda).*hcl3...
        -yoad.*J.*yyoil.*hcl4m3;
    
    dJ1=yoad.*hcl3;
    ddJ1dz=-dyadz2.*hcl3-yoad.*hcl4m3;
    ddJ1da=-dyada2.*hcl3-yoad.*yyoil.*hcl4m3;
    
    dJJ0=J.*hcl3;
    ddJJ0dz=dJdz.*hcl3-J.*hcl4m3;
    ddJJ0da=dJda.*hcl3-J.*yyoil.*hcl4m3;
    
%     dJJ2=(yy(:,idoil).^2-ya2.^2).*J./(h(:,idoil).^3);
    ddJJ2dz=(-2*ya2.*dyadz2.*J+yoad2.*dJdz).*hcl3...
        -yoad2.*J.*hcl4m3;
    ddJJ2da=(-2*ya2.*dyada2.*J+yoad2.*dJda).*hcl3...
        -yoad2.*J.*yyoil.*hcl4m3;
    
    dJ2=yoad2.*hcl3;
    ddJ2dz=-2*ya2.*dyadz2.*hcl3...
        -yoad2.*hcl4m3;
    ddJ2da=-2*ya2.*dyada2.*hcl3...
        -yoad2.*yyoil.*hcl4m3;
   
    J0=sum((dJ0(1:end-1,:)+dJ0(2:end,:))/2).*dyoil;
    dJ0dz=sum((ddJ0dz(1:end-1,:)+ddJ0dz(2:end,:))/2).*dyoil...
        +dybdz.*hb3-dyadz.*ha3;
    dJ0da=sum((ddJ0da(1:end-1,:)+ddJ0da(2:end,:))/2).*dyoil...
        +dybda.*hb3-dyada.*ha3;
    
%     JJ1=sum((dJJ1(1:end-1,:)+dJJ1(2:end,:))/2*dyoil2);
    dJJ1dz=sum((ddJJ1dz(1:end-1,:)+ddJJ1dz(2:end,:))/2).*dyoil...
        +(yb-ya).*Jb.*hb3.*dybdz;
    dJJ1da=sum((ddJJ1da(1:end-1,:)+ddJJ1da(2:end,:))/2).*dyoil...
        +(yb-ya).*Jb.*hb3.*dybda;
    
    J1=sum((dJ1(1:end-1,:)+dJ1(2:end,:))/2).*dyoil;
    dJ1dz=sum((ddJ1dz(1:end-1,:)+ddJ1dz(2:end,:))/2).*dyoil...
        +(yb-ya).*hb3.*dybdz;
    dJ1da=sum((ddJ1da(1:end-1,:)+ddJ1da(2:end,:))/2).*dyoil...
        +(yb-ya).*hb3.*dybda;
    
    JJ0=sum((dJJ0(1:end-1,:)+dJJ0(2:end,:))/2).*dyoil;
    dJJ0dz=sum((ddJJ0dz(1:end-1,:)+ddJJ0dz(2:end,:))/2).*dyoil...
        +Jb.*hb3.*dybdz;
    dJJ0da=sum((ddJJ0da(1:end-1,:)+ddJJ0da(2:end,:))/2).*dyoil...
        +Jb.*hb3.*dybda;
    
%     JJ2=sum((dJJ2(1:end-1,:)+dJJ2(2:end,:))/2*dyoil2);
    dJJ2dz=sum((ddJJ2dz(1:end-1,:)+ddJJ2dz(2:end,:))/2).*dyoil...
        +(yb.^2-ya.^2).*Jb.*hb3.*dybdz;
    dJJ2da=sum((ddJJ2da(1:end-1,:)+ddJJ2da(2:end,:))/2).*dyoil...
        +(yb.^2-ya.^2).*Jb.*hb3.*dybda;
    
    J2=sum((dJ2(1:end-1,:)+dJ2(2:end,:))/2).*dyoil;
    dJ2dz=sum((ddJ2dz(1:end-1,:)+ddJ2dz(2:end,:))/2).*dyoil...
        +(yb.^2-ya.^2).*hb3.*dybdz;
    dJ2da=sum((ddJ2da(1:end-1,:)+ddJ2da(2:end,:))/2).*dyoil...
        +(yb.^2-ya.^2).*hb3.*dybda;
    
    c1=(Pdoil-Pioil-12*mu_oilo.*JJ0)./J0;
    dc1dz=-12*mu_oilo./J0.*dJJ0dz-(Pdoil-Pioil-12*mu_oilo.*JJ0)./(J0.^2).*(dJ0dz);
    dc1da=-12*mu_oilo./J0.*dJJ0da-(Pdoil-Pioil-12*mu_oilo.*JJ0)./(J0.^2).*(dJ0da);
    
    dfgldz(ifoil)=Pioil.*dyadz-Pdoil.*dyadz-12*mu_oilo.*dJJ1dz-dJ1dz.*c1-dc1dz.*J1;
    dfgldalr(ifoil)=Pioil.*dyada-Pdoil.*dyada-12*mu_oilo.*dJJ1da-dJ1da.*c1-dc1da.*J1;
    
    dmgldz(ifoil)=ya.*dyadz.*Pioil-ya.*dyadz.*Pdoil-6*mu_oilo.*dJJ2dz-c1/2.*dJ2dz-1/2*dc1dz.*J2;
    dmgldalr(ifoil)=ya.*dyada.*Pioil-ya.*dyada.*Pdoil-6*mu_oilo.*dJJ2da-c1/2.*dJ2da-1/2*dc1da.*J2;
    
    dfgldpd(ifoil)=y2-ya-J1./J0;
    dfgldpi(ifoil)=ya-y1+J1./J0;
    dmgldpd(ifoil)=1/2*(y2^2-ya.^2)-1/2*J2./J0;
    dmgldpi(ifoil)=1/2*(ya.^2-y1^2)+1/2*J2./J0;
end
