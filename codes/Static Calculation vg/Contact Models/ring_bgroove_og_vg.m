function [fgl,mgl]=ring_bgroove_og_vg(Pi,Pd,gcl,dz_g,dz_go,zg,zr,zro,alr,alro,beta,betao,beta_th,y1,y2,hoil,mu_oil,dt)
Npe=length(zr);
Nby=50;

fgl=zeros(1,Npe);
mgl=zeros(1,Npe);

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

ifgas=minh>=hoil;
ngas=nnz(ifgas);
%gas region
if ngas
    ifgas1=ifgas&(abs(al)>1e-8);
    ifgas2=ifgas&(abs(al)<=1e-8);
    ngas1=nnz(ifgas1);
    ngas2=nnz(ifgas2);
    if ngas1
        h0g1=h0(ifgas1)-hoil(ifgas1);
        alg1=al(ifgas1);
        hc1=h0g1+alg1*y1;
        hc1=hc1.*(hc1>=1e-12)+1e-12*(hc1<1e-12);
        hc2=h0g1+alg1*y2;
        hc2=hc2.*(hc2>=1e-12)+1e-12*(hc2<1e-12);
        Pig1=Pi(ifgas1);
        Pdg1=Pd(ifgas1);
        
        h12p=(hc1+hc2);
        h12d=(hc1-hc2);
        fcom=h12d./alg1./h12p;
        mcom=2*(hc1.*hc2).^2.*log(hc1./hc2);
        mcomd=1./(alg1.^2.*h12p.*h12d);
        fgl(ifgas1)=-Pig1.*hc1.*fcom-Pdg1.*hc2.*fcom;
        mgl(ifgas1)=1/2*Pig1.*(mcom+hc1.*h12d.*(2*h0g1.*h12d-hc1.*h12p)).*mcomd...
            +1/2*Pdg1.*(-mcom+hc2.*h12d.*(2*h0g1.*h12d+hc2.*h12p)).*mcomd;
    end
    if ngas2
        Pig2=Pi(ifgas2);
        Pdg2=Pd(ifgas2);
        fgl(ifgas2)=1/2*(Pig2+Pdg2)*(y2-y1);
        mgl(ifgas2)=1/6*Pig2*(y2-y1)*(2*y1+y2)+1/6*Pdg2*(y2-y1)*(y1+2*y2);
    end
end

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
    
    yo=(hoil(ifoil)-h0oil)./(aloil);
    yo=min(yo,y2-1e-6);
    yo=max(yo,y1+1e-6);
    yb=y2*(aloil<0)+yo.*(aloil>0);
    %yb=((yb-y1)>1e-8).*yb+(y1+1e-8).*((yb-y1)<1e-8);    %yb cannot smaller than y1
    ya=y1*(aloil>0)+yo.*(aloil<0);    
    %ya=((y2-ya)>1e-8).*ya+(y2-1e-8).*((y2-ya)<1e-8);  %ya cannot excess y2
    ya2=repmat(ya,Nby+1,1);
    
    dyoil=(yb-ya)/Nby;
    dyoil2=repmat(dyoil,Nby,1);
    yyoil=repmat(ya,Nby+1,1)+repmat([0:Nby]',1,noil).*repmat(dyoil,Nby+1,1);
    h0_2oil=repmat(h0oil,Nby+1,1);
    al_2oil=repmat(aloil,Nby+1,1);
    hcl=h0_2oil+al_2oil.*yyoil;
    hcl=hcl.*(hcl>1e-7)+1e-7*(hcl<1e-7);
    hcl3=hcl.^(-3);
    yoad=yyoil-ya2;
    yoad2=yyoil.^2-ya2.^2;
    
    h0_2oilo=repmat(h0ooil,Nby+1,1);
    al_2oilo=repmat(alooil,Nby+1,1);
    hclo=h0_2oilo+al_2oilo.*yyoil;
    hclo=hclo.*(hclo>1e-7)+1e-7*(hclo<1e-7);
    
    dhdt=(hcl-hclo)/dt;
    dhdt=dhdt.*(dhdt<0);
    
    J=cumsum((dhdt(1:end-1,:)+dhdt(2:end,:))/2).*dyoil2;
    J=[zeros(1,noil);J];
    dJ0=hcl3;
    dJJ0=J.*hcl3;
    dJJ1=yoad.*J.*hcl3;
    dJ1=yoad.*hcl3;
    dJJ2=yoad2.*J.*hcl3;
    dJ2=yoad2.*hcl3;
    
    J0=sum((dJ0(1:end-1,:)+dJ0(2:end,:))/2).*dyoil;
    JJ1=sum((dJJ1(1:end-1,:)+dJJ1(2:end,:))/2).*dyoil;
    J1=sum((dJ1(1:end-1,:)+dJ1(2:end,:))/2).*dyoil;
    JJ0=sum((dJJ0(1:end-1,:)+dJJ0(2:end,:))/2).*dyoil;
    JJ2=sum((dJJ2(1:end-1,:)+dJJ2(2:end,:))/2).*dyoil;
    J2=sum((dJ2(1:end-1,:)+dJ2(2:end,:))/2).*dyoil;
    
    c1=(Pdoil-Pioil-12*mu_oilo.*JJ0)./J0;
    fgl(ifoil)=Pioil.*(ya-y1)+Pdoil.*(y2-ya)-12*mu_oilo.*JJ1-J1.*c1;
    mgl(ifoil)=1/2*(ya.^2-y1^2).*Pioil+1/2*(y2^2-ya.^2).*Pdoil-6*mu_oilo.*JJ2-c1/2.*J2;
end
